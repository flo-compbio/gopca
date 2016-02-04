# Copyright (c) 2015 Florian Wagner
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

import sys
import os
import logging
import re
import cPickle as pickle
import time
import hashlib
from copy import deepcopy
from collections import OrderedDict, Iterable
import datetime

import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import pearsonr

import goparser
from goparser import GOParser

from genometools.basic import GeneSetDB
from genometools.expression import ExpMatrix, ExpGene, ExpGenome

import gopca
from gopca import util
from gopca import enrichment
from gopca.enrichment import GSEAnalysis, GSEResult
from gopca import GOPCAConfig, GOPCASignature, GOPCAResult, GOPCARun

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
        GO-PCA configuration data.
    """

    def __init__(self, config):
        # store configuration
        assert isinstance(config, GOPCAConfig)
        self.config = deepcopy(config)

    @staticmethod
    def print_signatures(signatures, maxlength = 50, debug = False):
        """Print a list of signatures, sorted by their enrichment score.
        """
        a = sorted(range(len(signatures)),
                key = lambda i: -signatures[i].escore)

        for i in a:
            sig = signatures[i]
            sig_label = sig.get_label(max_name_length = maxlength,
                    include_pval = True)
            if debug:
                logger.debug(sig_label)
            else:
                logger.info(sig_label)

    @staticmethod
    def _read_expression(config):
        """Read expression data and perform variance filtering."""

        E = ExpMatrix.read_tsv(config.expression_file)
        # ExpMatrix constructor automatically sorts genes alphabetically

        p,n = E.shape
        logger.info('Expression matrix size: ' +
                '(p = %d genes) x (n = %d samples).', p,n)

        sel_var_genes = config.sel_var_genes
        if sel_var_genes > 0:
            # perform variance filtering

            if sel_var_genes >= p:
                # the number of genes to retain is equal or larger the number
                # of genes in the expresison matrix
                logger.warning('Variance filtering has no effect ' +
                        '(sel_var_genes = %d, p = %d).', sel_var_genes, p)
                return E

            #genes,samples,E = (exp.genes,exp.samples,exp.E)

            # select most variable genes
            var = np.var(E.X, axis=1, ddof=1)
            total_var = np.sum(var) # total sum of variance
            a = np.argsort(var)
            a = a[::-1]
            sel = np.zeros(p, dtype = np.bool_)
            sel[a[:sel_var_genes]] = True
            sel = np.nonzero(sel)[0]

            # report some information about the excluded genes
            lost_p = p - sel.size
            lost_var = total_var - np.sum(var[sel])
            logger.info('Selected the %d most variable genes ' +
                    '(excluded %.1f%% of genes, representing %.1f%% ' +
                    'of total variance).',
                    sel_var_genes, 100 * (lost_p / float(p)),
                    100 * (lost_var / total_var))

            # filtering
            E.genes = [E.genes[i] for i in sel]
            E.X = E.X[sel,:]

            p,n = E.shape
            logger.info('Expression matrix size, after variance filtering: ' +
                    'p = %d genes x n = %d samples.', p, n)

        return E

    @staticmethod
    def _read_gene_ontology(config):
        """Read the Gene Ontology data."""

        p_logger = logging.getLogger(goparser.__name__)
        p_logger.setLevel(logging.ERROR)
        P = GOParser()
        P.parse_ontology(config.gene_ontology_file,
                part_of_cc_only = config.go_part_of_cc_only)
        p_logger.setLevel(logging.NOTSET)

        return P

    @staticmethod
    def _estimate_n_components(config, X):
        """Estimate the number of non-trivial PCs using a permutation test."""

        assert isinstance(config, GOPCAConfig)
        assert isinstance(X, np.ndarray) 
        logger.info('Estimating the number of principal components ' +
                '(seed = %d)...', config.pc_seed)
        logger.debug('(permutations = %d, z-score threshold = %.1f)...',
                config.pc_permutations, config.pc_zscore_thresh)

        # perform PCA
        p,n = X.shape
        d_max = min(p, n - 1)
        M_pca = PCA(n_components = d_max)
        M_pca.fit(X.T)

        d = M_pca.explained_variance_ratio_
        logger.debug('Largest explained variance: %.2f', d[0])

        thresh = util.get_pc_explained_variance_threshold(X,
                config.pc_zscore_thresh, config.pc_permutations,
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
        a = sorted(range(q), key = lambda i: -enriched[i].escore)
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
        for i,g in enumerate(ranked_genes):
            if g not in genes_used:
                new_ranked_genes.append(g)
            elif i < L: # gene was already used, adjust L if necessary
                new_L -= 1
        ranked_genes = new_ranked_genes
        L = new_L

        ### start filtering

        # suppress logging messages from the enrichment module
        enr_logger = logging.getLogger(enrichment.__name__)
        enr_logger.setLevel(logging.ERROR)

        K_max = max([enr.K for enr in todo])
        p = len(ranked_genes)
        # initialize matrix for XL-mHG test
        mat = np.zeros((K_max + 1, p + 1), dtype = np.longdouble)
        while todo:
            most_enriched = todo[0]
            gs_id = most_enriched.gene_set.id

            # test if GO term is still enriched after removing all previously
            # used genes
            enr = M_enrich.get_enriched_gene_sets(ranked_genes, config.pval_thresh,
                    config.mHG_X_frac, config.mHG_X_min, L,
                    escore_pval_thresh = config.escore_pval_thresh,
                    gene_set_ids = [gs_id], mat = mat)
            assert len(enr) in [0,1]
            # enr will be an empty list if GO term does not meet the p-value
            # threshold

            todo = todo[1:] # remove the current gene set from the to-do list
            if not enr:
                continue
            elif config.escore_thresh is not None and \
                    enr[0].escore < config.escore_thresh:
                continue

            enr = enr[0]
            #print enr,'%d @ %d, s=%.1e' %(enr.k_n,enr.mHG_n,enr.stat)

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

        logger.info('Local filter: Kept %d / %d enriched gene sets.', \
                len(kept), len(enriched))

        return kept

    @staticmethod
    def _generate_signature(config, genes, X, pc, result):
        """Algorithm for generating a signature based on an enriched gene set.
        """

        # select genes above cutoff giving rise to XL-mHG test statistic
        enr_genes = result.genes[:result.k_n]

        # calculate average expression
        indices = np.int64([genes.index(g) for g in enr_genes])
        X_enr = X[indices,:]
        mean = np.mean(util.get_standardized_matrix(X_enr), axis = 0)

        # calculate seed based on the X genes most strongly correlated with the average
        corr = np.float64([pearsonr(mean, x)[0] for x in X_enr])
        a = np.argsort(corr)
        a = a[::-1]
        seed = np.mean(util.get_standardized_matrix(X_enr[a[:result.X],:]), 0)

        # select all other genes with correlation of at least sig_corr_thresh
        additional_indices = np.int64([i for i in a[result.X:]
                if pearsonr(seed, X_enr[i,:])[0] >= config.sig_corr_thresh])
        sel = np.r_[a[:result.X], additional_indices]
        sig_genes = [enr_genes[i] for i in sel]
        sig_X = X_enr[sel,:]

        return GOPCASignature(sig_genes, sig_X, pc, result)

    @staticmethod
    def _generate_pc_signatures(config, genes, E, M, W, pc):
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
        a = np.argsort(W[:,pc_index])
        if pc > 0:
            a = a[::-1]
        ranked_genes = [M.genome[i].name for i in a]

        # - find enriched gene sets using the XL-mHG test
        # - get_enriched_gene_sets() also calculates the enrichment score,
        #   but does not use it for filtering
        logger.debug('config: %f %d %d' %(config.mHG_X_frac, config.mHG_X_min, config.mHG_L))
        enriched = M.get_enriched_gene_sets(ranked_genes, config.pval_thresh,
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

        #logger.debug('-'*70)
        #logger.debug('All enriched gene sets:')
        #for enr in enriched:
        #    logger.debug(enr.get_pretty_format())
        #logger.debug('-'*70)

        # filter enriched gene sets (local filter)
        if not config.no_local_filter:
            enriched = GOPCA._local_filter(config, M, enriched, ranked_genes)

        # generate signatures
        signatures = []
        q = len(enriched)
        for j, enr in enumerate(enriched):
            signatures.append(GOPCA._generate_signature(config, genes, E, pc, enr))
        logger.info('Generated %d signatures based on the enriched gene sets.',
                q)

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
            term = go_parser.terms[term_id] # get the GOTerm object
            novel = True
            for t in (set([term_id]) | term.ancestors | term.descendants):
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
    ### end static functions

    ### public functions
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

    def run(self):
        """Run GO-PCA.

        Parameters
        ----------
        None

        Returns
        -------
        `gopca.GOPCARun`
            The GO-PCA run.
        """
        t0 = time.time()
        # check the configuration
        if not self.config.check():
            # problems with the configuration
            return None

        # get the timestamp
        timestamp = str(datetime.datetime.utcnow())
        logger.info('Timestamp: %s', timestamp)

        # make a copy of the configuration
        config = deepcopy(self.config)

        if config.gene_ontology_file is None and (not config.no_global_filter):
            # no ontology file => disable global filter
            logger.warning('Disabling global filter, since no gene ontology ' +
                    'file was provided.')
            config.set_param('no_global_filter', True)

        # read expression
        logger.info('Reading expression data...')
        hashval = util.get_file_md5sum(config.expression_file)
        logger.info('Expression file hash: %s', hashval)
        if config.expression_file_hash is not None:
            if config.expression_file_hash == hashval:
                logging.info('MD5 hash of expression file matches specified value.')
            else:
                logging.error('MD5 hash of expression file does not match!')
                return None
        else:
            config.set_param('expression_file_hash', hashval)
        E = self._read_expression(config)

        # read ontology
        go_parser = None
        if config.gene_ontology_file is not None:
            logger.info('Reading gene ontology...')
            logger.debug('part_of_cc_only: %s', str(config.go_part_of_cc_only))
            hashval = util.get_file_md5sum(config.gene_ontology_file)
            logger.info('Gene ontology file hash: %s', hashval)
            if config.gene_ontology_file_hash is not None:
                if config.gene_ontology_file_hash == hashval:
                    logging.info('MD5 hash of gene ontology file matches specified value.')
                else:
                    logging.error('MD5 hash of gene ontology file does not match specified value!')
                    return None
            else:
                config.set_param('gene_ontology_file_hash', hashval)
            go_parser = self._read_gene_ontology(config)

        # read gene sets
        logger.info('Reading gene sets...')
        hashval = util.get_file_md5sum(config.gene_set_file)
        logger.info('Gene set file hash: %s', hashval)
        if config.gene_set_file_hash is not None:
            if config.gene_set_file_hash == hashval:
                logging.info('MD5 hash of gene set file matches specified value.')
            else:
                logging.error('MD5 hash of gene set file does not match specified value!')
                return None
        else:
            config.set_param('gene_set_file_hash', hashval)
        D = GeneSetDB.read_tsv(config.gene_set_file)

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
                logger.error('The estimated number of non-trivial ' +
                        'principal components is zero!')

        else:
            d_max = min(E.p, E.n - 1)
            if config.n_components > d_max:
                logger.error('The number of PCs to test was specified as ' +
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
        M_enrich = GSEAnalysis(genome, D)

        # run PCA
        logger.info('Performing PCA...')
        M_pca = PCA(n_components = config.n_components)
        Y = M_pca.fit_transform(E.X.T)

        # output cumulative fraction explained for each PC
        frac = M_pca.explained_variance_ratio_
        cum_frac = np.cumsum(frac)
        logger.info('Cumulative fraction of variance explained by the first ' +
                '%d PCs: %.1f%%', config.n_components, 100 * cum_frac[-1])

        # generate signatures
        W = M_pca.components_.T
        final_signatures = []
        p = E.p
        var_expl = 0.0
        res_var = None
        for pc in range(config.n_components):
            var_expl += frac[pc]
            logger.info('')
            logger.info('-'*70)
            logger.info('PC %d explains %.1f%% of the variance.',
                    pc + 1, 100 * frac[pc])
            logger.info('The new cumulative fraction of variance explained ' +
                    'is %.1f%%.', 100 * var_expl)

            signatures = []
            signatures_dsc = self._generate_pc_signatures(config, E.genes, E.X, M_enrich, W, pc+1)
            signatures_asc = self._generate_pc_signatures(config, E.genes, E.X, M_enrich, W, -pc-1)
            signatures = signatures_dsc + signatures_asc

            logger.info('# signatures: %d', len(signatures))
            before = len(signatures)

            if not config.no_global_filter:
                signatures = self._global_filter(config, signatures, final_signatures, go_parser)
                logger.info('Global filter: kept %d / %d signatures.',
                        len(signatures),before)
        
            self.print_signatures(signatures, debug = True)
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
        run = GOPCARun(gopca.__version__, self.config, timestamp, result, run_time)

        if config.output_file is not None:
            logger.info('Storing GO-PCA run in file "%s"...',
                    config.output_file)
            run.write_pickle(config.output_file)

        return run
