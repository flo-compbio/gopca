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
from collections import OrderedDict

import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import pearsonr

import goparser
from goparser import GOParser

from genometools.expression import ExpMatrix

from gopca import util
from gopca import go_enrichment
from gopca.go_enrichment import GOEnrichmentAnalysis
from gopca import GOPCAConfig, GOPCASignature, GOPCAOutput

logger = logging.getLogger(__name__)

class GOPCA(object):
    """Class for performing GO-PCA.

    This class implements the GO-PCA algorithm, except for the GO enrichment
    testing, which is implemented by the `go_enrichment.GOEnrichmentAnalyis`
    class. The input data is provided as a `gopca.GOPCAConfig` object during
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
        self.__config = deepcopy(config)

    def __getattr__(self,name):
        """Custom attribute lookup function.

        We use this function to simplify access to GO-PCA parameters, which are
        stored in a `GOPCAConfig` object. Each (unknown) attribute name
        starting with an underscore gets mapped to the `GOPCAConfig` attribute
        of the same name, with the underscore removed.

        For example, this allows us to write, ``self._expression_file`` instead
        of ``self.__config.expression_file``.
        """
        # map attributes to parameters
        if name[0] == '_' and self.has_param(name[1:]):
            return getattr(self.__config,name[1:])
        raise AttributeError

    ### static methods
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
    def _test_file_hash(path, hashval):
        h = util.get_file_md5sum(path)
        return h == hashval

    @staticmethod
    def _read_expression(config):
        """Read expression data and perform variance filtering."""

        logger.info('Reading expression data...')
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
        if config.gene_ontology_file is None:
            raise AttributeError('No ontology file provided!')

        logger.info('Reading ontology...')
        logger.debug('part_of_cc_only: %s', str(config.go_part_of_cc_only))

        p_logger = logging.getLogger(goparser.__name__)
        p_logger.setLevel(logging.ERROR)
        P = GOParser()
        P.parse_ontology(config.gene_ontology_file,
                part_of_cc_only = config.go_part_of_cc_only)
        p_logger.setLevel(logging.NOTSET)

        return P

    @staticmethod
    def _read_go_annotations(config):
        """Read the GO annotations."""
        logger.info('Reading GO annotations...')
        go_annotations = util.read_go_annotations(config.go_annotation_file)
        return go_annotations

    @staticmethod
    def _estimate_n_components(config, X):
        """Estimate the number of non-trivial PCs using a permutation test."""

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
        config.set_param('n_components',d_est)

    @staticmethod
    def _local_filter(config, M_enrich, enriched_terms, ranked_genes):
        """Apply GO-PCA's "local" filter.
        
        Returns enrichments that pass the filter.
        """

        if len(enriched_terms) <= 1:
            return enriched_terms

        # sort enriched terms by enrichment
        q = len(enriched_terms)
        a = sorted(range(q), key = lambda i: -enriched_terms[i].escore)
        todo = [enriched_terms[i] for i in a]

        # keep the most enriched term
        most_enriched = todo[0]
        kept_terms = [most_enriched]
        todo = todo[1:]

        # exclude genes annotated with the most enriched term
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

        # start filtering
        enr_logger = logging.getLogger(go_enrichment.__name__)
        enr_logger.setLevel(logging.ERROR)
        K_max = max([enr.K for enr in todo])
        p = len(ranked_genes)
        mat = np.zeros((K_max + 1, p + 1), dtype=np.longdouble)
        while todo:
            most_enriched = todo[0]
            term_id = most_enriched.term[0]

            # test if GO term is still enriched after removing all previously
            # used genes
            enr = M_enrich.get_enriched_terms(ranked_genes, config.pval_thresh,
                    config.mHG_X_frac, config.mHG_X_min, L,
                    escore_pval_thresh = config.escore_pval_thresh,
                    selected_term_ids = [term_id], mat = mat)
            assert len(enr) in [0,1]
            # enr will be an empty list if GO term does not meet the p-value
            # threshold
            if enr:
                enr = enr[0]
                #print enr,'%d @ %d, s=%.1e' %(enr.mHG_k_n,enr.mHG_n,enr.mHG_s)

                # test fold enrichment threshold (if specified)
                still_enriched = False
                if config.escore_thresh is None:
                    still_enriched = True
                elif enr.escore >= config.escore_thresh:
                    still_enriched = True

                if still_enriched:
                    # if GO term is still considered enriched, keep it!
                    kept_terms.append(most_enriched)

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

            # next!
            todo = todo[1:]

        enr_logger.setLevel(logging.NOTSET)
        logger.info('Local filter: Kept %d / %d enriched terms.', \
                len(kept_terms), len(enriched_terms))

        return kept_terms

    @staticmethod
    def _generate_signature(config, genes, X, pc, enr):
        """Algorithm for generating a signature based on an enriched GO term.
        """

        # select genes above cutoff giving rise to XL-mHG test statistic
        enr_genes = enr.genes[:enr.mHG_k_n]

        # calculate average expression
        indices = np.int64([genes.index(g) for g in enr_genes])
        X_enr = X[indices,:]
        mean = np.mean(util.get_standardized_matrix(X_enr), axis = 0)

        # calculate seed based on the X genes most strongly correlated with the average
        corr = np.float64([pearsonr(mean, x)[0] for x in X_enr])
        a = np.argsort(corr)
        a = a[::-1]
        seed = np.mean(util.get_standardized_matrix(X_enr[a[:enr.X],:]), 0)

        # select all other genes with correlation of at least sig_corr_thresh
        additional_indices = np.int64([i for i in a[enr.X:]
                if pearsonr(seed, X_enr[i,:])[0] >= config.sig_corr_thresh])
        sel = np.r_[a[:enr.X], additional_indices]
        sig_genes = [enr_genes[i] for i in sel]
        sig_X = X_enr[sel,:]

        return GOPCASignature(sig_genes, sig_X, pc, enr)

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
        pc_index = abs(pc)-1
        a = np.argsort(W[:,pc_index])
        if pc > 0:
            a = a[::-1]
        ranked_genes = [genes[i] for i in a]

        # - find enriched GO terms using the XL-mHG test
        # - get_enriched_terms() also calculates the enrichment score,
        #   but does not use it for filtering
        logger.debug('config: %f %d %d' %(config.mHG_X_frac, config.mHG_X_min, config.mHG_L))
        enriched_terms = M.get_enriched_terms(ranked_genes, config.pval_thresh,
                config.mHG_X_frac, config.mHG_X_min, config.mHG_L,
                config.escore_pval_thresh)
        if not enriched_terms:
            return []

        # filter enriched GO terms by strength of enrichment
        # (if threshold is provided)
        if config.escore_thresh is not None:
            q_before = len(enriched_terms)
            enriched_terms = [t for t in enriched_terms
                    if t.escore >= config.escore_thresh]
            q = len(enriched_terms)
            logger.info('Kept %d / %d enriched terms with E-score >= %.1f',
                    q, q_before, config.escore_thresh)

        #logger.debug('-'*70)
        #logger.debug('All enriched terms:')
        #for enr in enriched_terms:
        #    logger.debug(enr.get_pretty_format())
        #logger.debug('-'*70)

        # filter enriched GO terms (local filter)
        if not config.no_local_filter:
            enriched_terms = GOPCA._local_filter(config, M, enriched_terms, ranked_genes)

        # generate signatures
        signatures = []
        q = len(enriched_terms)
        for j,enr in enumerate(enriched_terms):
            signatures.append(GOPCA._generate_signature(config, genes, E, pc, enr))
        logger.info('Generated %d signatures based on the enriched GO terms.',
                q)

        return signatures

    @staticmethod
    def _global_filter(config, new_signatures, previous_signatures, go_parser):
        """Apply GO-PCA's "global" filter.
        """
        if len(previous_signatures) == 0:
            return new_signatures

        kept_signatures = []
        previous_terms = set([sig.term[0] for sig in previous_signatures])
        for sig in new_signatures:
            term_id = sig.term[0]
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
                kept_signatures.append(sig)
        return kept_signatures

    @staticmethod
    def _get_config_dict(config):
        return config.get_dict()
    ### end static functions

    ### public functions
    def has_param(self, name):
        return self.__config.has_param(name)

    def get_param(self, name):
        return self.__config.get_param(name)

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
        self.__config.set_param(name, value)

    def run(self):
        """Run GO-PCA.

        Parameters
        ----------
        None

        Returns
        -------
        `gopca.GOPCAOutput`
            The GO-PCA output.
        """
        t0 = time.time()

        # check the configuration
        if not self.__config.check():
            # problems with the configuration
            return 1

        # make a copy of the configuration
        config = deepcopy(self.__config)

        if config.gene_ontology_file is None and (not config.no_global_filter):
            # no ontology file => disable global filter
            logger.warning('Disabling global filter, since no gene ontology ' +
                    'file was provided.')
            config.set_param('no_global_filter', True)

        # read expression
        hashval = util.get_file_md5sum(config.expression_file)
        if config.expression_file_hash is not None and \
                config.expression_file_hash != hashval:
            logging.error('MD5 hash of expression file does not match!')
            return 1
        else:
            config.expression_file_hash = hashval
        logger.info('Expression file hash: %s', config.expression_file_hash)
        E = self._read_expression(config)

        # determine mHG_L, if 0 or None
        if config.mHG_L is None:
            # None = "default" value
            config.set_param('mHG_L', int(len(exp.genes)/8.0))
        elif config.mHG_L == 0:
            # 0 = "disabled"
            config.set_param('mHG_L', len(exp.genes))

        # read ontology
        if config.gene_ontology_file is not None:
            hashval = util.get_file_md5sum(config.gene_ontology_file)
            if config.gene_ontology_file_hash is not None and \
                    config.gene_ontology_file_hash != hashval:
                logging.error('MD5 hash of gene ontology file does not match!')
                return 1
            else:
                config.gene_ontology_file_hash = hashval
            logger.info('Gene ontology file hash: %s',
                    config.gene_ontology_file_hash)
            go_parser = self._read_gene_ontology(config)

        # read GO annotations
        hashval = util.get_file_md5sum(config.go_annotation_file)
        if config.go_annotation_file_hash is not None and \
                config.go_annotation_file_hash != hashval:
            logging.error('MD5 hash of GO annotation file does not match!')
            return 1
        else:
            config.go_annotation_file_hash = hashval
        logger.info('GO annotation file hash: %s',
                config.go_annotation_file_hash)
        go_annotations = self._read_go_annotations(config)


        if config.n_components == 0:
            # estimate the number of non-trivial PCs using a permutation test
            self._estimate_n_components(config, E.X)
            if config.n_components == 0:
                # no non-trivial PCs => abort
                logger.error('The estimated number of non-trivial ' +
                        'principal components is zero!')
                raise ValueError

        logger.debug('-'*70)
        logger.debug('GO-PCA will be run with these parameters:')
        for d in config.get_param_strings():
            logger.debug(d)
        logger.debug('-'*70)

        # create GOEnrichment object
        M_enrich = GOEnrichmentAnalysis(E.genes, go_annotations)

        # perform PCA
        logger.info('Performing PCA...')
        M_pca = PCA(n_components = config.n_components)
        Y = M_pca.fit_transform(E.X.T)

        # output cumulative fraction explained for each PC
        frac = M_pca.explained_variance_ratio_
        cum_frac = np.cumsum(frac)
        logger.info('Cumulative fraction of variance explained by the first %d PCs: %.1f%%', \
                config.n_components, 100 * cum_frac[-1])

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
                    pc+1, 100*frac[pc])
            logger.info('The new cumulative fraction of variance explained ' +
                    'is %.1f%%.', 100*var_expl)

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
        logger.info('GO-PCA generated %d signatures:', len(final_signatures))
        self.print_signatures(final_signatures)

        S = np.float64([util.get_signature_expression(E.genes, E.X, sig.genes)
                for sig in final_signatures])

        # include the input data in the output data
        output = GOPCAOutput(self.__config, config, E.genes, E.samples, W, Y,
                final_signatures, S)

        t1 = time.time()
        logger.info('Total GO-PCA runtime: %.2f s.', t1-t0)

        return output


