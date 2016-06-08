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
from genometools.expression import filter as exp_filter
from genometools import enrichment
from genometools.enrichment import GSEAnalysis
from genometools.ontology import GeneOntology

import gopca
from . import GOPCAConfig, GOPCASignature, GOPCASignatureMatrix, GOPCARun
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
    matrix: `genometools.expression.ExpMatrix`
        The expression matrix.
    gene_set_db: `genometools.basics.GeneSetDB`
        See :attr:`gene_sets` attribute.
    ontology: `genometools.ontology.GeneOntology`, optional
        See :attr:`ontology` attribute.

    Attributes
    ----------
    config: `GOPCAConfig`
        GO-PCA configuration data.
    gene_set_db: `genometools.basics.GeneSetDB`
        The gene sets.
    ontology: `genometools.ontology.GeneOntology` or None
        The Gene Ontology.
    """
    def __init__(self, config, matrix, gene_set_db, ontology=None):
        # store configuration
        assert isinstance(config, GOPCAConfig)
        assert isinstance(gene_set_db, GeneSetDB)
        assert isinstance(matrix, ExpMatrix)
        if ontology is not None:
            assert isinstance(ontology, GeneOntology)
        
        self.config = config
        self.matrix = matrix
        self.gene_set_db = gene_set_db
        self.ontology = ontology

    @staticmethod
    def print_signatures(signatures, maxlength=50, debug=False):
        """Print a list of signatures, sorted by PC and enrichment score.
        """
        sig_sorted = sorted(signatures, \
                            key=lambda sig: [abs(sig.pc), -sig.pc, sig.escore])
        for sig in sig_sorted:
            sig_label = sig.get_label(max_name_length=maxlength,
                                      include_pval=True)
            if debug:
                logger.debug(sig_label)
            else:
                logger.info(sig_label)

    @staticmethod
    def _estimate_n_components(config, X):
        """Estimate the number of non-trivial PCs using a permutation test.

        This function also takes the ``pc_max`` parameter into account."""

        assert isinstance(config, GOPCAConfig)
        assert isinstance(X, np.ndarray) and X.ndim == 2
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

        logger.info('The estimated number of non-trivial PCs is %d.', d_est)

        if 0 < config.pc_max < d_est:
            d_est = config.pc_max
            logger.info('Limiting the number of PCs to test to %d.', d_est)

        config.set_param('n_components', d_est)

    @staticmethod
    def _local_filter(config, gse_analysis, enriched, ranked_genes):
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
        genes_used = set(most_enriched.ind_genes)
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


        # initialize matrix for XL-mHG test
        K_max = max([enr.K for enr in todo])
        p = len(ranked_genes)
        table = np.empty((K_max+1, p+1), dtype=np.longdouble)
        while todo:
            most_enriched = todo[0]
            gs_id = most_enriched.gene_set.id

            # test if GO term is still enriched after removing all previously
            # used genes
            enr = gse_analysis.get_enriched_gene_sets(
                ranked_genes, config.pval_thresh,
                config.mHG_X_frac, config.mHG_X_min, L,
                escore_pval_thresh=config.escore_pval_thresh,
                gene_set_ids=[gs_id], table=table)
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
            genes_used.update(most_enriched.ind_genes)
            new_ranked_genes = []
            new_L = L
            for i, g in enumerate(ranked_genes):
                if g not in genes_used:
                    new_ranked_genes.append(g)
                elif i < L:  # gene was already used, adjust L if necessary
                    new_L -= 1
            ranked_genes = new_ranked_genes
            L = new_L

        # stop suppressing log messages from the enrichment module
        enr_logger.setLevel(logging.NOTSET)

        logger.info('Local filter: Kept %d / %d enriched gene sets.',
                    len(kept), len(enriched))
        return kept

    @staticmethod
    def _generate_signature(config, matrix, pc, gse_result):
        """
        Algorithm for generating a signature based on an enriched gene set.
        """

        # select genes above cutoff giving rise to XL-mHG test statistic
        enr_genes = gse_result.genes_above_cutoff

        # calculate average expression (using standardized data)
        enr_matrix = matrix.loc[enr_genes].copy()
        enr_expr = enr_matrix.standardize_genes().mean(axis=0).values

        # calculate seed based on the X genes most strongly correlated with
        # the average
        corr = np.float64([pearsonr(enr_expr, x)[0] for x in enr_matrix.X])
        a = np.argsort(corr)
        a = a[::-1]
        enr_matrix_std = enr_matrix.standardize_genes()
        seed_expr = enr_matrix_std.iloc[a[:gse_result.X]].mean(axis=0).values

        # select all other genes with correlation of at least sig_corr_thresh
        additional_indices = np.int64(
            [i for i in a[gse_result.X:]
             if pearsonr(seed_expr, enr_matrix.iloc[i].values)[0] >=
             config.sig_corr_thresh])
        sel = np.r_[a[:gse_result.X], additional_indices]
        sig_matrix = enr_matrix.iloc[sel].copy()

        return GOPCASignature(pc, gse_result, sig_matrix)

    @staticmethod
    def _generate_pc_signatures(config, matrix, gse_analysis, W, pc):
        """Generate signatures for a specific principal component and ordering.

        The absolute value  of ``pc`` determines the principal component (PC).
        Genes are then ranked by their loadings for this PC. Whether this
        ranking is in ascending or descending order is determined by the sign
        of ``pc``: If it has a  positive sign, then the ranking will be in
        descending order (most positive loading values first). If it has a
        negative sign, then the ranking will be in ascending order (most
        negative loading values first).
        """
        assert isinstance(config, GOPCAConfig)
        assert isinstance(matrix, ExpMatrix)
        assert isinstance(gse_analysis, GSEAnalysis)
        assert isinstance(W, np.ndarray) and W.ndim == 2
        assert isinstance(pc, int) and pc != 0

        # rank genes by their PC loadings
        pc_index = abs(pc)-1
        a = np.argsort(W[:, pc_index])
        if pc > 0:
            # for positive pc values, use descending order
            a = a[::-1]
        ranked_genes = [matrix.index[i] for i in a]

        # - find enriched gene sets using the XL-mHG test
        # - get_enriched_gene_sets() also calculates the enrichment score,
        #   but does not use it for filtering
        logger.debug('config: %f %d %d',
                     config.mHG_X_frac, config.mHG_X_min, config.mHG_L)
        enriched = gse_analysis.get_enriched_gene_sets(
            ranked_genes, config.pval_thresh,
            config.mHG_X_frac, config.mHG_X_min, config.mHG_L,
            config.escore_pval_thresh)
        if not enriched:
            # no gene sets were found to be enriched
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

        # apply local filter (if enabled)
        if not config.no_local_filter:
            enriched = GOPCA._local_filter(config, gse_analysis,
                                           enriched, ranked_genes)

        # generate signatures
        signatures = []
        q = len(enriched)
        for j, enr in enumerate(enriched):
            signatures.append(
                GOPCA._generate_signature(config, matrix, pc, enr))
        logger.info(
            'Generated %d signatures based on the enriched gene sets.', q)

        return signatures

    @staticmethod
    def _global_filter(config, new_signatures, previous_signatures,
                       ontology=None):
        """Apply GO-PCA's "global" filter.

        """
        if len(previous_signatures) == 0:
            return new_signatures

        kept = []
        previous_gene_sets = set([sig.gene_set.id for sig in
                                 previous_signatures])
        for sig in new_signatures:
            gs_id = sig.gene_set.id

            test_gene_sets = {gs_id, }
            if ontology is not None:
                term = ontology[gs_id]  # get the GOTerm object
                test_gene_sets |= (term.ancestors | term.descendants)

            overlap = test_gene_sets & previous_gene_sets
            if overlap:
                logger.debug('Gene set "%s" filtered out.', sig.gene_set.name)
            else:
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

    def run(self):
        """Run GO-PCA.

        Parameters
        ----------

        Returns
        -------
        `GOPCARun` or None
            The GO-PCA run, or "None" if the run failed.
        """

        t0 = time.time()
        # check the configuration
        if not self.config.check_params():
            # problems with the configuration
            logger.error('Invalid configuration settings. '
                         'Aborting GO-PCA run.')
            return None

        # get the timestamp
        timestamp = str(datetime.datetime.utcnow())
        logger.info('Timestamp: %s', timestamp)

        # print some information
        p, n = self.matrix.shape
        logger.info('Size of expression matrix: ' +
                    'p=%d genes x n=%d samples.', p, n)

        # make a copy of the configuration
        config = copy.deepcopy(self.config)


        # get hash values
        config_hash = self.config.hash
        expression_hash = self.matrix.hash
        gene_sets_hash = self.gene_set_db.hash
        if self.ontology is not None:
            ontology_hash = self.ontology.hash
        else:
            ontology_hash = None
        logger.info('User configuration hash: %s', config_hash)
        logger.info('Expression matrix hash: %s', expression_hash)
        logger.info('Gene set hash: %s', gene_sets_hash)
        if ontology_hash is not None:
            logger.info('Gene set hash: %s', gene_sets_hash)

        # perform variance filtering (this creates a copy of E)
        matrix = self.matrix
        if config.sel_var_genes != 0:
            matrix = exp_filter.filter_variance(matrix, config.sel_var_genes)
            #logger.info('Expression matrix size, after variance filtering: ' +
            #            'p = %d genes x n = %d samples.', p, n)

        # determine mHG_L, if -1 or 0
        if config.mHG_L == -1:
            # -1 = determine L automatically => p / 8
            config.set_param('mHG_L', int(matrix.p / 8.0))
        elif config.mHG_L == 0:
            # 0 = "disable" effect of L => set it to the number of genes
            config.set_param('mHG_L', matrix.p)

        if config.n_components == -1:
            # estimate the number of non-trivial PCs using a permutation test
            self._estimate_n_components(config, matrix.X)
            if config.n_components == 0:
                logger.error('The estimated number of non-trivial '
                             'principal components is zero. '
                             'Aborting GO-PCA run.')
                return None

        else:
            d_max = min(matrix.p, matrix.n - 1)
            if config.n_components > d_max:
                logger.error('The number of PCs to test was specified as '
                             '%d, but the data only has %d PCs. '
                             'Aborting GO-PCA run.',
                             config.n_components, d_max)
                return None

        if config.n_components == 0:
            logger.error('No principal components to test.'
                         'Aborting GO-PCA run.')
            return None

        logger.debug('-'*70)
        logger.debug('GO-PCA will be run with these parameters:')
        for d in config.param_strings:
            logger.debug(d)
        logger.debug('-'*70)

        # create GSEAnalysis object
        genome = ExpGenome.from_gene_names(matrix.genes.tolist())
        gse_analysis = GSEAnalysis(genome, self.gene_set_db)

        # run PCA
        logger.info('Performing PCA...')
        pca = PCA(n_components=config.n_components)
        Y = pca.fit_transform(matrix.X.T)

        # output fraction of variance explained for the PCs tested
        frac = pca.explained_variance_ratio_
        cum_frac = np.cumsum(frac)
        logger.info('Fraction of variance explained by the first '
                    '%d PCs: %.1f%%', config.n_components, 100*cum_frac[-1])

        # generate signatures
        W = pca.components_.T
        final_signatures = []
        var_expl = 0.0
        for c in range(config.n_components):
            var_expl += frac[c]
            logger.info('')
            logger.info('-'*70)
            logger.info('PC %d explains %.1f%% of the variance.',
                        c+1, 100*frac[c])
            logger.info('The new cumulative fraction of variance explained '
                        'is %.1f%%.', 100*var_expl)

            signatures_dsc = self._generate_pc_signatures(
                config, matrix, gse_analysis, W, c+1)
            signatures_asc = self._generate_pc_signatures(
                config, matrix, gse_analysis, W, -(c+1))
            signatures = signatures_dsc + signatures_asc
            logger.info('# signatures: %d', len(signatures))

            # apply global filter (if enabled)
            if not config.no_global_filter:
                before = len(signatures)
                signatures = self._global_filter(
                    config, signatures, final_signatures, self.ontology)
                logger.info('Global filter: kept %d / %d signatures.',
                            len(signatures), before)
        
            self.print_signatures(signatures, debug=True)
            final_signatures.extend(signatures)
            logger.info('Total no. of signatures generated so far: %d',
                        len(final_signatures))

        logger.info('')
        logger.info('='*70)
        logger.info('GO-PCA generated %d signatures.',
                    len(final_signatures))

        # sort signatures?
        self.print_signatures(final_signatures)

        sig_matrix = GOPCASignatureMatrix(final_signatures, matrix.samples)
        t1 = time.time()
        exec_time = t1 - t0
        logger.info('This GO-PCA run took %.2f s.', exec_time)
        run = GOPCARun(gopca.__version__, timestamp,
                       self.config, config,
                       expression_hash, gene_sets_hash, ontology_hash,
                       matrix.genes, matrix.samples, W, Y,
                       sig_matrix, exec_time)

        return run
