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

"""Module containing the GOPCA class.
"""

import sys
import os
import logging
import re
import cPickle as pickle
import time
from copy import deepcopy
from collections import OrderedDict

import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import pearsonr

from gopca import common
from gopca import GOPCAInput

from goparser import GOParser
from gopca.go_enrichment import GOEnrichmentAnalysis

class GOPCA(object):
    """Class for performing GO-PCA.

    Parameters
    ----------
    input_: GOPCAInput object
        GO-PCA input data.
    """

    def __init__(self,input_):

        # configure logger
        self._logger = logging.getLogger(__name__)

        # store input data
        assert isinstance(input_,GOPCAInput)
        self.__raw_input = deepcopy(input_)

        # this is the "working copy"
        self.__input = deepcopy(self.__raw_input)

    def __getattr__(self,name):
        # map attributes to parameters
        if name[0] == '_' and name[1:] in self.__input.param_names:
            return getattr(self.__input,name[1:])
        raise AttributeError

    # logging convenience functions
    def _debug(self,s,*args):
        self._logger.debug(s,*args)

    def _info(self,s,*args):
        self._logger.info(s,*args)

    def _warning(self,s,*args):
        self._logger.warning(s,*args)

    def _error(self,s,*args):
        self._logger.error(s,*args)

    def _read_expression_data(self):
        """Read expression data and perform variance filtering."""

        self._info('Reading expression data...')
        genes,samples,E = common.read_expression(self._expression_file)

        p,n = E.shape
        self._info('Expression matrix size: ' +
                'p = %d genes x n = %d samples.', p,n)

        sel_var_genes = self._sel_var_genes
        if sel_var_genes > 0:
            # perform variance filtering

            if sel_var_genes >= p:
                # the number of genes to retain is equal or larger the number
                # of genes in the expresison matrix
                self._warning('Variance filtering has no effect ' +
                        '(sel_var_genes = %d, p = %d).', sel_var_genes, p)
                return genes,samples,E

            # rank genes by their variance
            var = np.var(E,axis=1,ddof=1)
            total_var = np.sum(var) # total sum of variance
            a = np.argsort(var)
            a = a[::-1]

            # select the genes with the largest variance
            p,n = E.shape
            sel = np.zeros(p,dtype=np.bool_)
            sel[a[:sel_var_genes]] = True
            sel = np.nonzero(sel)[0]

            # report some information about the excluded genes
            lost_p = p - sel.size
            lost_var = total_var - np.sum(var[sel])
            self._info('Selected the %d most variable genes ' +
                    '(excluded %.1f%% of genes, representing %.1f%% ' +
                    'of total variance).',
                    sel_var_genes, 100*(lost_p/float(p)),
                    100*(lost_var/total_var))

            # filtering
            genes = [genes[i] for i in sel]
            E = E[sel,:]

            p,n = E.shape
            self._info('Expression matrix size, after variance filtering: ' +
                    'p = %d genes x n = %d samples.', p, n)

        return genes,samples,E

    def _estimate_n_components(self,E):

        assert isinstance(E,np.ndarray) 
        self._info('Estimating the number of principal components ' +
                '(seed = %d)...', self._seed)
        self._debug('(permutations = %d, z-score threshold = %.1f)...',
                self._pc_permutations, self._pc_zscore_thresh)

        # perform PCA
        p,n = E.shape
        d_max = min(p,n-1)
        M_pca = PCA(n_components = d_max)
        M_pca.fit(E.T)

        d = M_pca.explained_variance_ratio_
        self._debug('Largest explained variance: %.2f', d[0])

        thresh = common.get_pc_explained_variance_threshold(E,
                self._pc_zscore_thresh, self._pc_permutations, self._seed)
        self._debug('Explained variance threshold: %.2f', thresh)
        d_est = np.sum(d >= thresh)

        self._info('The estimated number of PCs is %d.', d_est)
        self.set_param('n_components',d_est)

    def print_signatures(self,signatures):
        a = None
        maxlength = 50
        a = sorted(range(len(signatures)),key=lambda i: -signatures[i].escore)

        for i in a:
            sig = signatures[i]
            self.message(sig.get_label(max_name_length=maxlength,include_pval=True))

    def read_ontology(self):
        go_parser = None
        if self._ontology_file is not None:
            self._info('Reading ontology...')
            go_parser = GOParser()
            go_parser.parse_ontology(self._ontology_file,
                    part_of_cc_only = self._go_part_of_cc_only)
        return go_parser

    def read_go_annotations(self):
        self._info('Reading GO annotations...')
        go_annotations = common.read_go_annotations(self._go_annotation_file)
        return go_annotations

        #n_assign = sum(len(v) for k,v in self.annotations.iteritems())
        #self.message('(%d annotations) done!', n_assign)

    def set_param(self,name,value):
        setattr(self.__input,name,value)

    def run(self):
        """Run GO-PCA."""

        t0 = time.time()

        # make sure the input data is valid
        self.__raw_input.validate()

        # calculate and report input hash
        self.__raw_input.calculate_hash()
        self._info('GO-PCA input hash: %s', self.__raw_input.hash)

        # make a copy of the input data before changing anything
        self.__input = deepcopy(self.__raw_input)

        if self._ontology_file is None and (not self._disable_global_filter):
            # no ontology file => disable global filter
            self._warning('Disabling global filter, since no ontology file was provided.')
            self.set_param('disable_global_filter',True)

        # read expression
        genes, samples, E = self._read_expression_data()

        # determine mHG_L, if 0 or None
        if self._mHG_L is None:
            # None = "default" value
            self.set_param('mHG_L',int(len(genes)/8.0))
        elif self._mHG_L == 0:
            # 0 = "disabled"
            self.set_param('mHG_L',len(genes))

        # read ontology
        go_parser = self.read_ontology()

        # read GO annotations
        go_annotations = self.read_go_annotations()

        if self._n_components == 0:
            # estimate the number of non-trivial PCs using a permutation test
            self._estimate_n_components(E)
            if self._n_components == 0:
                # no non-trivial PCs => abort
                self._logger.error('The estimated number of non-trivial ' +
                        'principal components is zero!')
                raise ValueError


        # create GOEnrichment object
        M_enrich = GOEnrichmentAnalysis(genes,go_annotations)

        # perform PCA
        self._info('Performing PCA...')
        M_pca = PCA(n_components = self._n_components)
        Y = M_pca.fit_transform(E.T)

        # output cumulative fraction explained for each PC
        frac = M_pca.explained_variance_ratio_
        cum_frac = np.cumsum(frac)
        self._info('Cumulative fraction of variance explained by the first %d PCs: %.1f%%', \
                self._n_components,100*cum_frac[-1])

        # generate signatures
        W = M_pca.components_.T
        final_signatures = []
        p = len(genes)
        var_expl = 0.0
        res_var = None
        for pc in range(self._n_components):

            var_expl += frac[pc]
            self._info('')
            self._info('-'*70)
            self._info('PC %d explains %.1f%% of the variance.',
                    pc+1, 100*frac[pc])
            self._info('The new cumulative fraction of variance explained ' +
                    'is %.1f%%.', 100*var_expl)

            signatures = []
            #signatures_dsc = self.get_pc_signatures(M_enrich,W,pc+1)
            #signatures_asc = self.get_pc_signatures(M_enrich,W,-pc-1)
            #signatures = signatures_dsc + signatures_asc

            #self.info_('# signatures: %d',len(signatures))
            #before = len(signatures)

            #if not self.config.disable_global_filter:
            #    signatures = self.global_filter(signatures,final_signatures,self.GO)
            #    self.message('Global filter: kept %d / %d signatures.', len(signatures),before)
        
            self.print_signatures(signatures)
            final_signatures.extend(signatures)
            self._info('Total no. of signatures so far: %d',
                    len(final_signatures))

            pc += 1

        self.info_('')
        self.info_('='*70)
        self.info_('GO-PCA generated %d signatures:', len(final_signatures))
        self.print_signatures(final_signatures)

        S = np.float64([common.get_signature_expression(genes,E,sig.genes) for sig in final_signatures])

        # it's important that we generate a copy of self.config for the result
        result = GOPCAResult(deepcopy(self.input),self.genes,self.samples,W,Y,final_signatures,S)

        t1 = time.time()
        self._info('Total GO-PCA runtime: %.2fs', t1-t0)

        return result

    def local_filter(self,M_enrich,enriched_terms,ranked_genes):
        # implements GO-PCA's "local" filter
        # returns enrichments that pass the filter

        if len(enriched_terms) <= 1:
            return enriched_terms

        # sort enriched terms by enrichment
        q = len(enriched_terms)
        a = sorted(range(q),key=lambda i:-enriched_terms[i].escore)
        todo = [enriched_terms[i] for i in a]

        # keep the most enriched term
        most_enriched = todo[0]
        kept_terms = [most_enriched]
        todo = todo[1:]

        # exclude genes annotated with the most enriched term
        genes_used = set(most_enriched.genes)
        new_ranked_genes = []
        L = self.L
        new_L = L
        for i,g in enumerate(ranked_genes):
            if g not in genes_used:
                new_ranked_genes.append(g)
            elif i < L: # gene was already used, adjust L if necessary
                new_L -= 1
        ranked_genes = new_ranked_genes
        L = new_L

        # start filtering
        K_max = max([enr.K for enr in todo])
        p = len(ranked_genes)
        mat = np.zeros((K_max+1,p+1),dtype=np.longdouble)
        while todo:
            most_enriched = todo[0]
            term_id = most_enriched.term[0]

            # test if GO term is still enriched after removing all previously used genes
            enr = M_enrich.get_enriched_terms(ranked_genes,self.pval_thresh,\
                    self.X_frac,self.X_min,L,\
                    escore_pval_thresh=self.escore_pval_thresh,selected_term_ids=[term_id],\
                    mat=mat,quiet=True)
            assert len(enr) in [0,1]
            if enr: # enr will be an empty list if GO term does not meet the p-value threshold
                enr = enr[0]
                #print enr,'%d @ %d, s=%.1e' %(enr.mHG_k_n,enr.mHG_n,enr.mHG_s)

                # test fold enrichment threshold (if specified)
                still_enriched = False
                if self.escore_thresh is None:
                    still_enriched = True
                elif enr.escore >= self.escore_thresh:
                    still_enriched = True

                if still_enriched:
                    # if GO term is still considered enriched, keep it!
                    kept_terms.append(most_enriched)
                    # next, exclude selected genes from further analysis: 1) adjust L 2) update set of excluded genes
                    genes_used.update(most_enriched.genes) # add all genes of the GO term to set of used genes
                    new_ranked_genes = []
                    new_L = L
                    for i,g in enumerate(ranked_genes):
                        if g not in genes_used:
                            new_ranked_genes.append(g)
                        elif i < L: # gene was already used, adjust L if necessary
                            new_L -= 1
                    ranked_genes = new_ranked_genes
                    L = new_L

            # next!
            todo = todo[1:]

        self.message('Filtering: Kept %d / %d enriched terms.', \
                len(kept_terms),len(enriched_terms))
        return kept_terms

    def global_filter(self,new_signatures,previous_signatures,GO):
        if len(previous_signatures) == 0:
            return new_signatures
        kept_signatures = []
        previous_terms = set([sig.term[0] for sig in previous_signatures])
        for sig in new_signatures:
            term_id = sig.term[0]
            term = GO.terms[term_id]
            novel = True
            for t in set([term_id]) | term.ancestors | term.descendants:
                if t in previous_terms:
                    novel = False
                    break
            if novel:
                kept_signatures.append(sig)
        return kept_signatures

    def get_pc_signatures(self,M,W,pc):
        """
        Generate GO-PCA signatures for a specific PC and a specific ranking of loadings (ascending or descending).
        The absolute value of the 'pc' parameter determines the principal component. Genes are then ranked by their loadings for this PC.
        Whether this ranking is in ascending or descending order is determined by the sign of the 'pc' parameter.
        -> If the 'pc' parameter has a positive sign, then the ranking will be in descending order (most positive loading values first)
        -> If the 'pc' parameter has a negative sign, then the ranking will be in ascending order (most negative loading values first)
        """

        # rank genes by their PC loadings
        pc_index = abs(pc)-1
        a = np.argsort(W[:,pc_index])
        if pc > 0:
            a = a[::-1]
        ranked_genes = [self.genes[i] for i in a]

        # find enriched GO terms using the XL-mHG test
        # get_enriched_terms also calculates the enrichment score, but does not use it for filtering
        enriched_terms = M.get_enriched_terms(ranked_genes,self.pval_thresh,self.X_frac,self.X_min,self.L,self.escore_pval_thresh)
        if not enriched_terms:
            return []

        # filter enriched GO terms by strength of enrichment (if threshold is provided)
        if self.escore_thresh is not None:
            q_before = len(enriched_terms)
            enriched_terms = [t for t in enriched_terms if t.escore >= self.escore_thresh]
            q = len(enriched_terms)
            self.message('Enrichment filter: Kept %d / %d enriched terms with enrichment score >= %.1fx.', \
                    q,q_before,self.escore_thresh)

        # filter enriched GO terms (local filter)
        if not self.config.disable_local_filter:
            q_before = len(enriched_terms)
            enriched_terms = self.local_filter(M,enriched_terms,ranked_genes)
            q = len(enriched_terms)
            self.message('Local filter: Kept %d / %d enriched terms.', q,q_before)

        # generate signatures
        signatures = []
        q = len(enriched_terms)
        for j,enr in enumerate(enriched_terms):
            signatures.append(self.get_signature(pc,enr))
        self.message('Generated %d GO-PCA signatures based on the enriched GO terms.', q)

        return signatures

    def get_signature(self,pc,enr):
        """ Generate signature seed by selecting the X genes most correlated with average expression. """

        # select genes above cutoff giving rise to XL-mHG test statistic
        enr_genes = enr.genes[:enr.mHG_k_n]

        # calculate average expression
        indices = np.int64([self.genes.index(g) for g in enr_genes])
        E_enr = self.E[indices,:]
        mean = np.mean(common.get_standardized_matrix(E_enr),axis=0)

        # calculate seed based on the X genes most strongly correlated with the average
        corr = np.float64([pearsonr(mean,e)[0] for e in E_enr])
        a = np.argsort(corr)
        a = a[::-1]
        seed = np.mean(common.get_standardized_matrix(E_enr[a[:enr.X],:]),0)

        # select all other genes with correlation of at least sig_corr_thresh
        additional_indices = np.int64([i for i in a[enr.X:] if pearsonr(seed,E_enr[i,:])[0] > self.sig_corr_thresh])
        sel = np.r_[a[:enr.X],additional_indices]
        sig_genes = [enr_genes[i] for i in sel]
        sig_E = E_enr[sel,:]

        return GOPCASignature(sig_genes,sig_E,pc,enr)



