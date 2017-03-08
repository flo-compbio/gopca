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
from collections import Iterable, OrderedDict
from pkg_resources import parse_version

import numpy as np
import sklearn
from sklearn.decomposition import PCA
from scipy.stats import pearsonr

from genometools.basic import GeneSetCollection
from genometools.expression import ExpProfile, ExpMatrix, ExpGene, ExpGenome
from genometools import enrichment
from genometools.enrichment import RankBasedGSEResult, \
                                   GeneSetEnrichmentAnalysis
from genometools.ontology import GeneOntology

import gopca
from . import GOPCAParams, GOPCAConfig, \
              GOPCASignature, GOPCASignatureMatrix, GOPCARun
from . import util

logger = logging.getLogger(__name__)


class GOPCA(object):
    """Class for performing GO-PCA.

    This class implements the GO-PCA algorithm. (The GO enrichment testing
    is implemented in the `enrichment.GeneSetEnrichmentAnalysis` class of
    the `genometools` package). The input data consists of an expression
    matrix (`genometools.expression.ExpMatrix`) and a list of GO-PCA
    "configurations" (`GOPCAConfig`), i.e., pairs of parameter settings and
    gene set collections.

    Parameters
    ----------
    matrix : `genometools.expression.ExpMatrix`
        See :attr:`matrix` attribute.
    configs : Iterable of `GOPCAConfig`
        See :attr:`configs` attribute.
    num_components : int, optional
        See :attr:`num_components` attribute. [0]
    pc_seed : int, optional
        See :attr:`pc_seed` attribute. [0]
    pc_num_permutations : int, optional
        See :attr:`pc_num_permutations` attribute. [15]
    pc_zscore_thresh : float, optional
        See :attr:`pc_zscore_thresh` attribute. [2.0]
    pc_max_components : int, optional
        See :attr:`pc_max_components` attribute. [0]
    verbose : bool, optional
        See :attr:`verbose` attribute. [False]

    Attributes
    ----------
    matrix : `genometools.expression.ExpMatrix`
        The expression matrix.
    configs : list of `GOPCAConfig`
        The list of GO-PCA configurations. Each configuration consists of
        gene sets (represented by a `GOPCAGeneSets` instance) along with a set
        of GO-PCA parameters (`GOPCAParams`) to use for testing those gene
        sets.
    num_components : int
        The number of principal components to test. If set 0, the number is
         determined automatically using a permutation-based algorithm.
    pc_seed : int
        The random number generator seed, used to generate the permutations
        for automatically determining the number of principal components to
        test.
    pc_num_permutations : int
        The number of permutations to used for automatically determining the
        number of principal components to test.
    pc_zscore_thresh : float
        The z-score threshold used for automatically determining the number of
        principal components (PC) to test. First, the fraction of variance
        explained by the first PC in each permuted dataset is calculated.
        Then, the mean and standard deviation of those values are used to
        calculate a z-score for the fraction of variance explained by each PC
        in the real dataset. All PCs with a z-score above the specified
        threshold are tested.
    pc_max_components : int
        The maximum number of principal components (PCs) to test (only relevant
        when the algorithm for automatically determining the number of PCs
        to test is used. For testing a fixed number of PCs, set the
        :attr:`num_components` attribute to a non-zero value.
    verbose : bool
        If set to ``True``, generate more verbose output.
    """
    __param_defaults = OrderedDict([
        ('num_components', 0),  # 0 = automatic
        ('pc_seed', 0),
        ('pc_num_permutations', 15),
        ('pc_zscore_thresh', 2.0),
        ('pc_max_components', 0),  # 0 = no maximum
    ])
    """Global GO-PCA parameter default values."""

    @staticmethod
    def get_param_defaults():
        return GOPCA.__param_defaults.copy()

    def __init__(self, matrix, configs, **kwargs):

        assert isinstance(matrix, ExpMatrix)
        assert isinstance(configs, Iterable)

        num_components = kwargs.pop(
            'num_components', self.__param_defaults['num_components'])

        pc_seed = kwargs.pop(
            'pc_seed', self.__param_defaults['num_components'])

        pc_num_permutations = kwargs.pop(
            'pc_num_permutations',
            self.__param_defaults['pc_num_permutations'])

        pc_zscore_thresh = kwargs.pop(
            'pc_zscore_thresh', self.__param_defaults['pc_zscore_thresh'])

        pc_max_components = kwargs.pop(
            'pc_max_components', self.__param_defaults['pc_max_components'])

        verbose = kwargs.pop('verbose', False)

        assert isinstance(num_components, (int, np.integer))
        assert isinstance(pc_seed, (int, np.integer))
        assert isinstance(pc_num_permutations, (int, np.integer))
        assert isinstance(pc_zscore_thresh, (float, np.float))
        assert isinstance(pc_max_components, (int, np.integer))
        assert isinstance(verbose, bool)

        self.matrix = matrix
        self.configs = list(configs)

        self.num_components = int(num_components)
        self.pc_seed = int(pc_seed)
        self.pc_num_permutations = int(pc_num_permutations)
        self.pc_zscore_thresh = float(pc_zscore_thresh)
        self.pc_max_components = int(pc_max_components)

        self.verbose = verbose

        # make sure configs have the right type
        for conf in self.configs:
            assert isinstance(conf, GOPCAConfig)

    @classmethod
    def simple_setup(cls, matrix, params, gene_sets, gene_ontology=None,
                     **kwargs):
        """Initialize GO-PCA instance with only one collection of gene sets.

        """
        # TODO: finish docstring
        assert isinstance(matrix, ExpMatrix)
        assert isinstance(params, GOPCAParams)
        assert isinstance(gene_sets, GeneSetCollection)
        if gene_ontology is not None:
            assert isinstance(gene_ontology, GeneOntology)

        configs = [GOPCAConfig(params, gene_sets, gene_ontology)]
        return cls(matrix, configs, **kwargs)

    def __repr__(self):
        return '<%s instance (hash="%s")>' % \
               (self.__class__.__name__, self.hash)


    def __str__(self):
        return '<%s instance (hash="%s")>' % \
               (self.__class__.__name__, self.hash)

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return self.__dict__ == other.__dict__
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def hash(self):
        data_str = ';'.join(
            repr(v) for v in [self.configs, self.matrix])
        data = data_str.encode('UTF-8')
        return str(hashlib.md5(data).hexdigest())

    @property
    def X(self):
        return self.matrix.X

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
    def get_pc_explained_variance_threshold(X, z, t, seed):
        # TODO: Write docstring.

        # RandomizedPCA does not work in Scikit-learn 0.14.1,
        # but it works in Scikit-learn 0.16.1
        use_old_randomized_pca = False
        sklearn_pv = parse_version(sklearn.__version__)

        use_old_pca = False 
        if sklearn_pv < parse_version('0.18.0'):
            use_old_pca = True
            if sklearn_pv >= parse_version('0.16.1'):
                # old randomized PCA implementation 
                logger.debug('Using old scikit-learn randomized PCA implementation.')
                from sklearn.decomposition import RandomizedPCA as PCA
            else:
                # no randomized PCA implementation available
                logger.debug('No randomized PCA implementation available.')
                from sklearn.decomposition import PCA
        else:
            # use new randomized PCA implementation
            from sklearn.decomposition import PCA
        
        # initialize random number generator
        np.random.seed(seed)

        # do permutations
        p, n = X.shape
        d_max_null = np.empty(t, dtype=np.float64)
        X_perm = np.empty((p, n), dtype=np.float64)
        if use_old_pca:
            M_null = PCA(n_components=1)
        else:
            M_null = PCA(n_components=1, svd_solver='randomized')
        
        for j in range(t):
            for i in range(p):
                X_perm[i, :] = X[i, np.random.permutation(n)]

            M_null.fit(X_perm.T)
            d_max_null[j] = M_null.explained_variance_ratio_[0]

        # calculate z-score threshold
        mean_null = np.mean(d_max_null)
        std_null = np.std(d_max_null, ddof=1)
        thresh = mean_null + z * std_null

        return thresh

    def estimate_num_components(self):
        """Estimate the number of non-trivial PCs using a permutation test.

        """
        # TODO: finish docstring
        logger.info('Estimating the number of principal components '
                    '(seed = %d)...', self.pc_seed)
        logger.debug('(permutations = %d, z-score threshold = %.1f)...',
                     self.pc_num_permutations, self.pc_zscore_thresh)

        # perform PCA
        p, n = self.matrix.shape
        d_max = min(p, n-1)
        M_pca = PCA(n_components=d_max)
        M_pca.fit(self.matrix.X.T)

        d = M_pca.explained_variance_ratio_
        logger.debug('Largest explained variance: %.2f', d[0])

        thresh = self.get_pc_explained_variance_threshold(
            self.X, self.pc_zscore_thresh, self.pc_num_permutations,
            self.pc_seed)
        logger.debug('Explained variance threshold: %.2f', thresh)
        d_est = np.sum(d >= thresh)

        logger.info('The estimated number of non-trivial PCs is %d.', d_est)

        return d_est

    @staticmethod
    def _local_filter(params, gse_analysis, enriched, ranked_genes,
                      verbose=False):
        """Apply GO-PCA's "local" filter.
        
        Returns the enriched gene sets that passed the filter.
        """

        assert isinstance(params, GOPCAParams)
        assert isinstance(gse_analysis, GeneSetEnrichmentAnalysis)
        assert isinstance(enriched, Iterable)
        assert isinstance(ranked_genes, Iterable)
        assert isinstance(verbose, bool)

        msg = logger.debug
        if verbose:
            msg = logger.info

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
        L = params.mHG_L
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
            enr = gse_analysis.get_rank_based_enrichment(
                ranked_genes, params.pval_thresh,
                params.mHG_X_frac, params.mHG_X_min, L,
                adjust_pval_thresh=False,
                escore_pval_thresh=params.escore_pval_thresh,
                gene_set_ids=[gs_id], table=table,
                exact_pval='if_necessary')
            assert len(enr) in [0, 1]
            # enr will be an empty list if GO term does not meet the p-value
            # threshold

            todo = todo[1:]  # remove the current gene set from the to-do list
            if not enr:
                continue
            elif params.escore_thresh is not None and \
                    enr[0].escore < params.escore_thresh:
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

        return kept

    @staticmethod
    def _generate_signature(matrix, params, pc, gse_result,
                            standardize=False, verbose=False):
        """Generate a signature based on an enriched gene set.
        """
        assert isinstance(matrix, ExpMatrix)
        assert isinstance(params, GOPCAParams)
        assert isinstance(pc, int)
        assert isinstance(gse_result, RankBasedGSEResult)
        assert isinstance(standardize, bool)
        assert isinstance(verbose, bool)

        # select genes above cutoff giving rise to XL-mHG test statistic
        enr_genes = gse_result.genes_above_cutoff

        # calculate average expression
        enr_matrix = matrix.loc[enr_genes].copy()

        if standardize:
            enr_matrix.standardize_genes(inplace=True)

        # use the average expression of all genes above the XL-mHG cutoff as
        # a "seed"
        seed = ExpProfile(enr_matrix.mean(axis=0))

        # rank all genes by their correlation with the seed, and select only
        # those with correlation ">=" params.sig_corr_thresh, but no fewer
        # than params.min_sig_genes

        # calculate seed based on the X genes most strongly correlated with
        # the average
        corr = np.float64([pearsonr(seed.values, x)[0] for x in enr_matrix.X])
        a = np.argsort(corr)
        a = a[::-1]

        # determine the number of genes to include
        num_genes = max(np.sum(corr >= params.sig_corr_thresh),
                        params.sig_min_genes)

        sig_matrix = enr_matrix.iloc[a[:num_genes]].copy()

        return GOPCASignature(pc, gse_result, seed, sig_matrix)

    @staticmethod
    def _generate_pc_signatures(matrix, params, gse_analysis, W, pc,
                                standardize=False, verbose=False):
        """Generate signatures for a specific principal component and ordering.

        The absolute value  of ``pc`` determines the principal component (PC).
        Genes are then ranked by their loadings for this PC. Whether this
        ranking is in ascending or descending order is determined by the sign
        of ``pc``: If it has a  positive sign, then the ranking will be in
        descending order (most positive loading values first). If it has a
        negative sign, then the ranking will be in ascending order (most
        negative loading values first).
        """
        assert isinstance(matrix, ExpMatrix)
        assert isinstance(params, GOPCAParams)
        assert isinstance(gse_analysis, GeneSetEnrichmentAnalysis)
        assert isinstance(W, np.ndarray) and W.ndim == 2
        assert isinstance(pc, int) and pc != 0
        assert isinstance(standardize, bool)
        assert isinstance(verbose, bool)

        msg = logger.debug
        if verbose:
            msg = logger.info

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

        # suppress logging messages from genometools.enrichment module
        enr_logger = logging.getLogger(enrichment.__name__)
        enr_logger.setLevel(logging.ERROR)

        logger.debug('config: %f %d %d',
                     params.mHG_X_frac, params.mHG_X_min, params.mHG_L)
        enriched = gse_analysis.get_rank_based_enrichment(
            ranked_genes, params.pval_thresh,
            params.mHG_X_frac, params.mHG_X_min, params.mHG_L,
            adjust_pval_thresh=False,
            escore_pval_thresh=params.escore_pval_thresh,
            exact_pval='if_significant')
        if not enriched:
            # no gene sets were found to be enriched
            return []

        # stop suppressing logging messages from genometools.enrichment module
        enr_logger.setLevel(logging.NOTSET)

        # filter enriched GO terms by strength of enrichment
        # (if threshold is provided)
        if params.escore_thresh is not None:
            q_before = len(enriched)
            enriched = [enr for enr in enriched
                        if enr.escore >= params.escore_thresh]
            q = len(enriched)
            msg('Kept %d / %d enriched gene sets with E-score >= %.1f',
                q, q_before, params.escore_thresh)

        # apply local filter (if enabled)
        if not params.no_local_filter:
            q_before = len(enriched)
            enriched = GOPCA._local_filter(params, gse_analysis,
                                           enriched, ranked_genes)
            q = len(enriched)
            msg('Local filter: Kept %d / %d enriched gene sets.', q, q_before)

        # generate signatures
        signatures = []
        q = len(enriched)
        for j, enr in enumerate(enriched):
            signatures.append(
                GOPCA._generate_signature(
                    matrix, params, pc, enr,
                    standardize=standardize, verbose=verbose))
        msg('Generated %d signatures based on the enriched gene sets.', q)

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
        """Perform GO-PCA.

        Parameters
        ----------

        Returns
        -------
        `GOPCARun` or None
            The GO-PCA run, or ``None`` if the run failed.
        """
        t0 = time.time()  # remember the start time
        timestamp = str(datetime.datetime.utcnow())  # timestamp for the run


        ### Phase 1: Make sure all configurations are valid
        all_configs_valid = True
        for config in self.configs:
            if not config.user_params.check_params():
                # problems with the configuration
                all_configs_valid = False
            config.finalize_params(self.matrix.p)
            if not config.params.check_params():
                all_configs_valid = False

        if not all_configs_valid:
            logger.error('Invalid configuration settings. '
                         'Aborting GO-PCA run.')
            return None

        # print some information
        p, n = self.matrix.shape
        logger.info('Timestamp: %s', timestamp)
        logger.info('Size of expression matrix: ' +
                    'p=%d genes x n=%d samples.', p, n)

        # Report hash values for expression matrix and configurations
        expression_hash = self.matrix.hash
        logger.info('Expression matrix hash: %s', expression_hash)
        config_hashes = []
        for i, config in enumerate(self.configs):
            config_hashes.append(config.hash)
            logger.info('Configuration #%d hash: %s', i+1, config_hashes[-1])


        ### Phase 2: Determine the number of principal components
        num_components = self.num_components
        if num_components == 0:
            # estimate the number of non-trivial PCs using a permutation test
            num_components = self.estimate_num_components()
            if num_components == 0:
                logger.error('The estimated number of non-trivial '
                             'principal components is 0. '
                             'Aborting GO-PCA run.')
                return None
            if 0 < self.pc_max_components < num_components:
                num_components = self.pc_max_components
                logger.info('Limiting the number of PCs to test to %d.', num_components)


        else:
            # determine the total number of principal components
            # (i.e., the number of dimensions spanned by the data)
            max_components = min(self.matrix.p, self.matrix.n - 1)
            if self.num_components > max_components:
                logger.error('The number of PCs to test was specified as '
                             '%d, but the data spans only %d dimensions. '
                             'Aborting GO-PCA run.',
                             num_components, max_components)
                return None

        if num_components == 0:
            logger.error('No principal components to test.'
                         'Aborting GO-PCA run.')
            return None


        ### Phase 3: Perform PCA
        logger.info('Performing PCA...')
        pca = PCA(n_components=num_components)
        Y = pca.fit_transform(self.matrix.X.T)

        # output fraction of variance explained for the PCs tested
        frac = pca.explained_variance_ratio_
        cum_frac = np.cumsum(frac)
        logger.info('Fraction of total variance explained by the first '
                    '%d PCs: %.1f%%', num_components, 100 * cum_frac[-1])


        ### Phase 4: Run GO-PCA for each configuration supplied
        enr_logger = logging.getLogger(enrichment.__name__)

        genome = ExpGenome.from_gene_names(self.matrix.genes.tolist())
        W = pca.components_.T  # the loadings matrix

        msg = logger.debug
        if self.verbose:
            # enable more verbose "INFO" messages
            msg = logger.info

        all_signatures = []
        for k, config in enumerate(self.configs):

            logger.info('Generating GO-PCA signatures for configuration '
                        '%d...', k+1)

            # create GeneSetEnrichmentAnalysis object
            enr_logger.setLevel(logging.ERROR)
            gse_analysis = GeneSetEnrichmentAnalysis(genome, config.gene_sets)
            enr_logger.setLevel(logging.NOTSET)

            # generate signatures
            final_signatures = []
            var_expl = 0.0
            for d in range(num_components):
                var_expl += frac[d]
                msg('')
                msg('-'*70)
                msg('PC %d explains %.1f%% of the variance.',
                    d+1, 100*frac[d])
                msg('The new cumulative fraction of variance explained '
                    'is %.1f%%.', 100*var_expl)

                signatures_dsc = self._generate_pc_signatures(
                    self.matrix, config.params, gse_analysis, W, d+1)
                signatures_asc = self._generate_pc_signatures(
                    self.matrix, config.params, gse_analysis, W, -(d+1))
                signatures = signatures_dsc + signatures_asc
                msg('# signatures: %d', len(signatures))

                # apply global filter (if enabled)
                if not config.params.no_global_filter:
                    before = len(signatures)
                    signatures = self._global_filter(
                        config.params, signatures, final_signatures,
                        config.gene_ontology)
                    msg('Global filter: kept %d / %d signatures.',
                        len(signatures), before)

                # self.print_signatures(signatures, debug=True)
                final_signatures.extend(signatures)
                msg('Total no. of signatures generated so far: %d',
                    len(final_signatures))

            logger.info('')
            logger.info('='*70)
            logger.info('GO-PCA for configuration #%d generated %d '
                        'signatures.', k+1, len(final_signatures))
            logger.info('-'*70)
            self.print_signatures(final_signatures)
            logger.info('='*70)
            logger.info('')
            all_signatures.extend(final_signatures)


        ### Phase 5: Generate signature matrix and return a `GOPCARun` instance
        sig_matrix = GOPCASignatureMatrix.from_signatures(all_signatures)
        t1 = time.time()
        exec_time = t1 - t0
        logger.info('This GO-PCA run took %.2f s.', exec_time)
        gopca_run = GOPCARun(sig_matrix,
                             gopca.__version__, timestamp, exec_time,
                             expression_hash, config_hashes,
                             self.matrix.genes, self.matrix.samples, W, Y)

        return gopca_run
