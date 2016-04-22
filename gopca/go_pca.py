# Copyright (c) 2015, 2016 Florian Wagner
#
# This file is part of GO-PCA.
#
# GO-PCA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Module containing the `GOPCA` class.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

# import sys
# import os
import logging
# import re
# import cPickle as pickle
import time
import hashlib
import copy
import datetime

import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import pearsonr

from goparser import GOParser

from genometools.basic import GeneSetDB
from genometools.expression import ExpMatrix, ExpGene, ExpGenome
from genometools import enrichment
from genometools.enrichment import GSEAnalysis

import gopca
from . import GOPCAConfig, GOPCASignature, GOPCAResult, GOPCARun
from . import util

logger = logging.getLogger(__name__)


class GOPCA(object):
    """Class for performing GO-PCA.

    This class implements the GO-PCA algorithm, except for the GO enrichment
    testing, which is implemented by the `enrichment.GSEAnalysis` class.
    The input data is provided as a `gopca.GOPCAConfig` object during
    instantiation, and GO-PCA is run from start to finish by calling the `run`
    method.

    Parameters
    ----------
    config: `go_pca.GOPCAConfig`
        See :attr:`config` attribute.
    gene_sets: `genometools.basics.GeneSetDB`
        See :attr:`gene_sets` attribute.
    go_parser: `goparser.GOParser`, optional
        See :attr:`go_parser` attribute.

    Attributes
    ----------
    config: `GOPCAConfig`
        GO-PCA configuration data.
    gene_sets: `genometools.basics.GeneSetDB`
        The gene sets.
    go_parser: `goparser.GOParser` or None
        The gene ontology.
    """

    def __init__(self, config, gene_sets, go_parser=None):
        # store configuration
        assert isinstance(config, GOPCAConfig)
        assert isinstance(gene_sets, GeneSetDB)
        if go_parser is not None:
            assert isinstance(go_parser, GOParser)
        
        self.config = config
        self.gene_sets = gene_sets
        self.go_parser = go_parser

    @staticmethod
    def print_signatures(signatures, maxlength=50, debug=False):
        """Print a list of signatures, sorted by their enrichment score.
        """
        a = sorted(range(len(signatures)),
                   key=lambda i: -signatures[i].escore)

        for i in a:
            sig = signatures[i]
            sig_label = sig.get_label(max_name_length=maxlength,
                                      include_pval=True)
            if debug:
                logger.debug(sig_label)
            else:
                logger.info(sig_label)

    @staticmethod
    def _variance_filter(config, E):
        """Perform variance filtering. Return a copy of E."""
        assert isinstance(config, GOPCAConfig)
        assert isinstance(E, ExpMatrix)

        sel_var_genes = config.sel_var_genes
        if sel_var_genes == 0:
            # variance filter is disabled
            return E

        p, n = E.shape

        if sel_var_genes > p:
            logger.warning('Variance filter parameter G = %d has no effect, '
                           'since there are only p = %d genes.'
                           % (sel_var_genes, p))
            return E.copy()

        # select most variable genes
        var = np.var(E.X, axis=1, ddof=1)
        total_var = np.sum(var)  # total sum of variance
        a = np.argsort(var)
        a = a[::-1]
        sel = np.zeros(p, dtype=np.bool_)
        sel[a[:sel_var_genes]] = True
        sel = np.nonzero(sel)[0]

        # report some information about the excluded genes
        lost_p = p - sel.size
        lost_var = total_var - np.sum(var[sel])
        logger.info('Selected the %d most variable genes '
                    '(excluded %.1f%% of genes, representing %.1f%% '
                    'of total variance).',
                    sel_var_genes, 100 * (lost_p / float(p)),
                    100 * (lost_var / total_var))

        # filtering (this creates a copy)
        E = E.iloc[sel]

        p, n = E.shape
        logger.info('Expression matrix size, after variance filtering: ' +
                    'p = %d genes x n = %d samples.', p, n)

        return E

    @staticmethod
    def _estimate_n_components(config, X):
        """Estimate the number of non-trivial PCs using a permutation test."""

        assert isinstance(config, GOPCAConfig)
        assert isinstance(X, np.ndarray) 
        logger.info('Estimating the number of principal components '
                    '(seed = %d)...', config.pc_seed)
        logger.debug('(permutations = %d, z-score threshold = %.1f)...',
                     config.pc_permutations, config.pc_zscore_thresh)

        # perform PCA
        p, n = X.shape
        d_max = min(p, n-1)
        M_pca = PCA(n_components=d_max)
        M_pca.fit(X.T)

        d = M_pca.explained_variance_ratio_
        logger.debug('Largest explained variance: %.2f', d[0])

        thresh = util.get_pc_explained_variance_threshold(
            X, config.pc_zscore_thresh, config.pc_permutations,
            config.pc_seed)
        logger.debug('Explained variance threshold: %.2f', thresh)
        d_est = np.sum(d >= thresh)

        logger.info('The estimated number of PCs is %d.', d_est)

        if config.pc_max > 0 and d_est > config.pc_max:
            d_est = config.pc_max
            logger.info('Limiting the number of PCs to test to %d.', d_est)

        config.set_param('n_components', d_est)

    @staticmethod
    def _local_filter(config, M_enrich, enriched, ranked_genes):
        """Apply GO-PCA's "local" filter.
        
        Returns the enriched gene sets that passed the filter.
        """
        if len(enriched) <= 1:
            return enriched

        # sort enriched gene sets by E-score (in descending order)
        q = len(enriched)
        a = sorted(range(q), key=lambda i: -enriched[i].escore)
        todo = [enriched[i] for i in a]

        # keep the most enriched gene set
        most_enriched = todo[0]
        kept = [most_enriched]
        todo = todo[1:]

        # exclude all genes contained in the most enriched gene set
        genes_used = set(most_enriched.genes)
        new_ranked_genes = []
        L = config.mHG_L
        new_L = L
        for i, g in enumerate(ranked_genes):
            if g not in genes_used:
                new_ranked_genes.append(g)
            elif i < L:  # gene was already used, adjust L if necessary
                new_L -= 1
        ranked_genes = new_ranked_genes
        L = new_L

        # start filtering

        # suppress logging messages from the enrichment module
        enr_logger = logging.getLogger(enrichment.__name__)
        enr_logger.setLevel(logging.ERROR)

        K_max = max([enr.K for enr in todo])
        p = len(ranked_genes)
        # initialize matrix for XL-mHG test
        mat = np.empty((K_max+1, p+1), dtype=np.longdouble)
        while todo:
            most_enriched = todo[0]
            gs_id = most_enriched.gene_set.id

            # test if GO term is still enriched after removing all previously
            # used genes
            enr = M_enrich.get_enriched_gene_sets(
                ranked_genes, config.pval_thresh,
                config.mHG_X_frac, config.mHG_X_min, L,
                escore_pval_thresh=config.escore_pval_thresh,
                gene_set_ids=[gs_id], mat=mat)
            assert len(enr) in [0, 1]
            # enr will be an empty list if GO term does not meet the p-value
            # threshold

            todo = todo[1:]  # remove the current gene set from the to-do list
            if not enr:
                continue
            elif config.escore_thresh is not None and \
                    enr[0].escore < config.escore_thresh:
                continue

            # enr = enr[0]
            # print enr,'%d @ %d, s=%.1e' %(enr.k_n,enr.mHG_n,enr.stat)

            # keep the gene set
            kept.append(most_enriched)

            # next, exclude selected genes from further analysis:
            # 1) update set of used (excluded) genes 2) adjust L
            genes_used.update(most_enriched.genes)
            new_ranked_genes = []
            new_L = L
            for i, g in enumerate(ranked_genes):
                if g not in genes_used:
                    new_ranked_genes.append(g)
                elif i < L: # gene was already used, adjust L if necessary
                    new_L -= 1
            ranked_genes = new_ranked_genes
            L = new_L

        # stop suppressing log messages from the enrichment module
        enr_logger.setLevel(logging.NOTSET)

        logger.info('Local filter: Kept %d / %d enriched gene sets.',
                    len(kept), len(enriched))

        return kept

    @staticmethod
    def _generate_signature(config, genes, samples, X, pc, result):
        """
        Algorithm for generating a signature based on an enriched gene set.
        """

        # select genes above cutoff giving rise to XL-mHG test statistic
        enr_genes = result.genes[:result.k_n]

        # calculate average expression
        indices = np.int64([genes.index(g) for g in enr_genes])
        X_enr = X[indices, :]
        mean = np.mean(util.get_standardized_matrix(X_enr), axis=0)

        # calculate seed based on the X genes most strongly correlated with
        # the average
        corr = np.float64([pearsonr(mean, x)[0] for x in X_enr])
        a = np.argsort(corr)
        a = a[::-1]
        seed = np.mean(util.get_standardized_matrix(X_enr[a[:result.X], :]), 0)

        # select all other genes with correlation of at least sig_corr_thresh
        additional_indices = np.int64(
            [i for i in a[result.X:]
             if pearsonr(seed, X_enr[i, :])[0] >= config.sig_corr_thresh])
        sel = np.r_[a[:result.X], additional_indices]
        sig_genes = [enr_genes[i] for i in sel]
        sig_X = X_enr[sel, :]

        return GOPCASignature(sig_genes, samples, sig_X, pc, result)

    @staticmethod
    def _generate_pc_signatures(config, genes, samples, E, M, W, pc):
        """Generate signatures for a specific principal component and ordering.

        The absolute value of ``pc`` determines the principal component (PC).
        Genes are then ranked by their loadings for this PC. Whether this
        ranking is in ascending or descending order is determined by the sign
        of ``pc``: If it has a  positive sign, then the ranking will be in
        descending order (most positive loading values first). If it has a
        negative sign, then the ranking will be in ascending order (most
        negative loading values first).
        """

        # rank genes by their PC loadings
        pc_index = abs(pc) - 1
        a = np.argsort(W[:, pc_index])
        if pc > 0:
            a = a[::-1]
        ranked_genes = [M.genome[i].name for i in a]

        # - find enriched gene sets using the XL-mHG test
        # - get_enriched_gene_sets() also calculates the enrichment score,
        #   but does not use it for filtering
        logger.debug('config: %f %d %d'
                     % (config.mHG_X_frac, config.mHG_X_min, config.mHG_L))
        enriched = M.get_enriched_gene_sets(
            ranked_genes, config.pval_thresh,
            config.mHG_X_frac, config.mHG_X_min, config.mHG_L,
            config.escore_pval_thresh)
        if not enriched:
            return []

        # filter enriched GO terms by strength of enrichment
        # (if threshold is provided)
        if config.escore_thresh is not None:
            q_before = len(enriched)
            enriched = [enr for enr in enriched
                        if enr.escore >= config.escore_thresh]
            q = len(enriched)
            logger.info('Kept %d / %d enriched gene sets with E-score >= %.1f',
                        q, q_before, config.escore_thresh)

        # logger.debug('-'*70)
        # logger.debug('All enriched gene sets:')
        # for enr in enriched:
        #    logger.debug(enr.get_pretty_format())
        # logger.debug('-'*70)

        # filter enriched gene sets (local filter)
        if not config.no_local_filter:
            enriched = GOPCA._local_filter(config, M, enriched, ranked_genes)

        # generate signatures
        signatures = []
        q = len(enriched)
        for j, enr in enumerate(enriched):
            signatures.append(
                GOPCA._generate_signature(config, genes, samples, E, pc, enr))
        logger.info(
            'Generated %d signatures based on the enriched gene sets.', q)

        return signatures

    @staticmethod
    def _global_filter(config, new_signatures, previous_signatures, go_parser):
        """Apply GO-PCA's "global" filter.

        Can only be applied for GO-derived gene sets.
        """
        if len(previous_signatures) == 0:
            return new_signatures

        kept = []
        previous_terms = set([sig.gene_set.id for sig in previous_signatures])
        for sig in new_signatures:
            term_id = sig.gene_set.id
            term = go_parser.terms[term_id]  # get the GOTerm object
            novel = True
            for t in ({term_id, } | term.ancestors | term.descendants):
                if t in previous_terms:
                    logger.debug('GO term "%s" filtered out due to "%s".',
                                 go_parser.terms[term_id].name,
                                 go_parser.terms[t].name)
                    novel = False
                    break
            if novel:
                kept.append(sig)
        return kept

    @staticmethod
    def _get_config_dict(config):
        return config.get_dict()
    # end static functions

    # public functions
    def has_param(self, name):
        return self.config.has_param(name)

    def get_param(self, name):
        return self.config.get_param(name)

    def set_param(self, name, value):
        """Set a GO-PCA parameter.

        Parameters
        ----------
        name: str
            The name of the parameter.
        value: ?
            The value of the parameter.

        Returns
        -------
        None
        """
        self.config.set_param(name, value)

    def run(self, E):
        """Run GO-PCA.

        Parameters
        ----------
        E: `genometools.expression.ExpMatrix`
            The expression matrix.

        Returns
        -------
        `gopca.GOPCARun` or None
            The GO-PCA run, or "None" if the run failed.
        """
        assert isinstance(E, ExpMatrix)

        t0 = time.time()
        # check the configuration
        if not self.config.check():
            # problems with the configuration
            return None

        # get the timestamp
        timestamp = str(datetime.datetime.utcnow())
        logger.info('Timestamp: %s', timestamp)

        # make a copy of the configuration
        config = copy.deepcopy(self.config)

        if self.go_parser is None and (not config.no_global_filter):
            # no ontology data => disable global filter
            logger.warning('Disabling global filter, since no gene ontology '
                           'data was provided.')
            config.set_param('no_global_filter', True)

        # generate expression hash
        # logger.info('Reading expression data...')
        # hashval = util.get_file_md5sum(config.expression_file)
        expression_hash = str(hashlib.md5(str(hash(E))).hexdigest())
        logger.info('Expression data hash value: %s', expression_hash)

        # generate ontology hash
        # => not implemented
        ontology_hash = None

        # generate gene sets hash
        gene_sets_hash = str(hashlib.md5(
            str(hash(self.gene_sets))).hexdigest())
        logger.info('Gene set data hash value: %s', gene_sets_hash)

        # perform variance filtering (this creates a copy of E)
        E = self._variance_filter(config, E)

        # determine mHG_L, if -1 or 0
        if config.mHG_L == -1:
            # -1 = determine L automatically => p / 8
            config.set_param('mHG_L', int(E.p / 8.0))
        elif config.mHG_L == 0:
            # 0 = "disable" effect of L => set it to the number of genes
            config.set_param('mHG_L', E.p)

        if config.n_components == -1:
            # estimate the number of non-trivial PCs using a permutation test
            self._estimate_n_components(config, E.X)
            if config.n_components == 0:
                logger.error('The estimated number of non-trivial '
                             'principal components is zero!')

        else:
            d_max = min(E.p, E.n - 1)
            if config.n_components > d_max:
                logger.error('The number of PCs to test was specified as '
                             '%d, but the data only has %d PCs.',
                             config.n_components, d_max)
                return None

        if config.n_components == 0:
            logger.error('No principal components to test.')
            return None

        logger.debug('-'*70)
        logger.debug('GO-PCA will be run with these parameters:')
        for d in config.get_param_strings():
            logger.debug(d)
        logger.debug('-'*70)

        # create GSEAnalysis object
        genome = ExpGenome([ExpGene(g) for g in E.genes])
        M_enrich = GSEAnalysis(genome, self.gene_sets)

        # run PCA
        logger.info('Performing PCA...')
        M_pca = PCA(n_components=config.n_components)
        Y = M_pca.fit_transform(E.X.T)

        # output cumulative fraction explained for each PC
        frac = M_pca.explained_variance_ratio_
        cum_frac = np.cumsum(frac)
        logger.info('Cumulative fraction of variance explained by the first '
                    '%d PCs: %.1f%%', config.n_components, 100*cum_frac[-1])

        # generate signatures
        W = M_pca.components_.T
        final_signatures = []
        # p = E.p
        var_expl = 0.0
        # res_var = None
        for pc in range(config.n_components):
            var_expl += frac[pc]
            logger.info('')
            logger.info('-'*70)
            logger.info('PC %d explains %.1f%% of the variance.',
                        pc + 1, 100 * frac[pc])
            logger.info('The new cumulative fraction of variance explained '
                        'is %.1f%%.', 100 * var_expl)

            signatures_dsc = self._generate_pc_signatures(
                config, E.genes, E.samples, E.X, M_enrich, W, pc+1)
            signatures_asc = self._generate_pc_signatures(
                config, E.genes, E.samples, E.X, M_enrich, W, -pc-1)
            signatures = signatures_dsc + signatures_asc

            logger.info('# signatures: %d', len(signatures))
            before = len(signatures)

            if not config.no_global_filter:
                signatures = self._global_filter(
                    config, signatures, final_signatures, self.go_parser)
                logger.info('Global filter: kept %d / %d signatures.',
                            len(signatures),before)
        
            self.print_signatures(signatures, debug=True)
            final_signatures.extend(signatures)
            logger.info('Total no. of signatures so far: %d',
                        len(final_signatures))

            pc += 1

        logger.info('')
        logger.info('='*70)
        logger.info('GO-PCA generated %d signatures.',
                    len(final_signatures))

        # sort signatures?
        self.print_signatures(final_signatures)

        S = np.float64([util.get_signature_expression(E.genes, E.X, sig.genes)
                        for sig in final_signatures])

        result = GOPCAResult(config, E.genes, E.samples, W, Y, final_signatures, S)
        t1 = time.time()
        run_time = t1 - t0
        logger.info('This GO-PCA run took %.2f s.', run_time)
        run_config = copy.deepcopy(self.config)
        run = GOPCARun(gopca.__version__, run_config,
                       expression_hash, gene_sets_hash, ontology_hash,
                       timestamp, run_time, result)

        return run
