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

import sys
import os
import csv
import logging
import gzip
import cPickle as pickle
from math import ceil

import numpy as np
from scipy.stats import hypergeom

from genometools import misc
import xlmhg

class mHGTermResult(object):
    """
    Stores mHG result for one particular term.
    """

    #Note: Change this so that it inherits from class mHGResult
    #     (only additional attributes: term, genes).

    def __init__(self,term,pval,genes,ranks,N,X,L,mHG_n=None,mHG_k_n=None,mHG_s=None):

        ranks = np.int32(ranks)
        ranks.flags.writeable=False # this makes ranks.data hashable

        self.term = term # 4-tuple (id,source,collection,name)
        self.genes = tuple(genes) # genes corresponding to the "1's"
        self.pval = pval # XL-mHG p-value
        self.ranks = ranks # ranks of the genes in the list
        self.N = N # total number of genes
        self.X = X # XL-mHG "X" parameter
        self.L = L # XL-mHG "L" parameter

        # data related to the XL-mHG test statistic
        self.mHG_n = mHG_n
        self.mHG_k_n = mHG_k_n
        self.mHG_s = mHG_s

        # strength of enrichment
        self.escore_pval_thresh = None
        self.escore = None

    def __repr__(self):
        return '<mHGTermResult: %s (pval=%.1e; X=%d; L=%d; N=%d); genes hash=%d; positions hash=%d)>' \
                %('/'.join(self.term),self.pval,self.X,self.L,self.N,hash(self.genes),hash(self.ranks.data))

    def __str__(self):
        return '<mHGTermResult of term "%s" (%s; %d genes) with p-value = %.1e (X=%d, L=%d, N=%d)>' \
                %(str(self.term[3]),self.term[0],len(self.genes),self.pval,self.X,self.L,self.N)

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self,other):
        if type(self) != type(other):
            return False
        elif repr(self) == repr(other):
            return True
        else:
            return False

    @property
    def k(self):
        return self.mHG_k_n

    @property
    def K(self):
        return self.ranks.size

    @property
    def n(self):
        return self.mHG_n

    @property
    def fold_enrichment(self):
        return self.k / (self.K * (self.n/float(self.N)))

    def calculate_escore(self,pval_thresh=1.0):
        """ Calculate XL-mHG enrichment score.  """
        N = self.N
        K = self.K
        X = self.X
        L = self.L
        ranks = self.ranks
        
        if K == 0 or L == N or K < X:
            return 0

        k_max = 0
        fe_max = 0.0
        k = 1
        pval = 1.0

        while k <= K and ranks[k-1] < L:
            if k >= X:
                n = ranks[k-1] + 1
                if pval_thresh == 1.0 or hypergeom.sf(k-1,N,K,n) <= pval_thresh:
                    fe = k / (K * (n / float(N)))
                    if fe >= fe_max:
                        fe_max = fe
                        k_max = k
            k += 1

        self.escore_pval_thresh = pval_thresh
        self.escore = fe_max

    def get_pretty_format(self,omit_param=True,max_name_length=0):
        term_name = self.term[3]
        if max_name_length > 0 and len(term_name) > max_name_length:
            assert max_name_length >= 3
            term_name = term_name[:(max_name_length-3)] + '...'
        term_str = term_name + ' (%d)' %(len(self.genes))
        param_str = ''
        if not omit_param:
            param_str = ' [X=%d,L=%d,N=%d]' %(self.X,self.L,self.N)
        escore_str = ''
        if self.escore is not None:
            escore_str = ', e=%.1fx' %(self.escore)
        details = ', p=%.1e%s%s' %(self.pval,escore_str,param_str)
        return '%s%s' %(term_str,details)
        
    def get_pretty_GO_format(self,GO,omit_acc=False,omit_param=True,max_name_length=0):
        # accepts a GOParser object ("GO")
        term = GO.terms[self.term[0]]
        term_name = term.get_pretty_format(omit_acc=omit_acc,max_name_length=max_name_length)
        term_str = term_name + ' (%d)' %(len(self.genes))
        param_str = ''
        if not omit_param:
            param_str = ' [X=%d,L=%d,N=%d]' %(self.X,self.L,self.N)
        escore_str = ''
        if self.escore is not None:
            escore_str = ', e=%.1fx' %(self.escore)
        details = ', p=%.1e%s%s' %(self.pval,escore_str,param_str)
        return '%s%s' %(term_str,details)

class GOEnrichment(object):

    def __init__(self,genes,annotations,logger):

        a = np.lexsort([genes])
        self.genes = [genes[i] for i in a]
        self.terms = sorted(annotations.keys(), key=lambda x:x[0])
        self.logger = logger.getChild('GOEnrich')
        #self.logger.propagate = False
        #term_ids = [t[0] for t in self.terms] # self.term_ids?
        #self.terms = [GO.terms[id_] for id_ in self.term_ids] # 4-tuples
        p = len(genes)
        m = len(annotations)
        self.A = np.zeros((p,m),dtype=np.uint8)
        for j,t in enumerate(self.terms):
            for g in annotations[t]:
                try:
                    idx = misc.bisect_index(self.genes,g)
                except ValueError:
                    pass
                else:
                    self.A[idx,j] = 1

    # logging convenience functions
    def message(self,s,*args):
        self.logger.info(s,*args)

    def warning(self,s,*args):
        self.logger.warning(s,*args)

    def error(self,s,*args):
        self.logger.error(s,*args)

    def get_enriched_terms(self,ranked_genes,pval_thresh,X_frac,X_min,L,escore_pval_thresh=None,selected_term_ids=[],mat=None,quiet=False,verbose=False):
        """
        Tests GO term enrichment of either all terms or the terms specified by ``selected_term_ids''.
        """

        log_level = logging.INFO
        if quiet:
            log_level = logging.WARNING
        elif verbose:
            log_level = logging.ERROR
        self.logger.setLevel(log_level)

        genes = self.genes
        terms = self.terms
        A = self.A

        if escore_pval_thresh is None:
            escore_pval_thresh = pval_thresh

        # test only some terms?
        if selected_term_ids:
            term_ids = [t[0] for t in terms]
            term_indices = np.int64([misc.bisect_index(term_ids,t) for t in selected_term_ids])
            terms = [terms[i] for i in term_indices]
            A = A[:,term_indices] # not a view!

        # sort rows in annotation matrix (and exclude genes not in the ranking)
        gene_indices = np.int64([misc.bisect_index(genes,g) for g in ranked_genes])
        A = A[gene_indices,:] # not a view either!

        # determine largest K
        K_lim = np.sum(A[:L,:],axis=0,dtype=np.int64)
        K_rem = np.sum(A[L:,:],axis=0,dtype=np.int64)
        K = K_lim + K_rem
        K_max = np.amax(K)

        # prepare matrix for XL-mHG p-value calculation
        p,m = A.shape
        if mat is None:
            mat = np.empty((K_max+1,p+1),dtype=np.longdouble)

        # find enriched GO terms
        self.message('Testing %d terms for enrichment...', m)
        self.message('(N = %d, X_frac = %.2f, X_min = %d, L = %d; K_max = %d)', \
                    len(ranked_genes),X_frac,X_min,L,K_max)

        enriched_terms = []
        tested = 0
        tested_mHG = 0
        for j in range(m):
            v = np.ascontiguousarray(A[:,j]) # copy

            # determine significance of enrichment using XL-mHG test
            # (only if there are at least X_min genes annotated with this term before the L'th threshold)
            #if K_lim[j] >= X_min:
            if K[j] >= X_min:
                tested += 1
                # determine term-specific X (based on K[j])
                X = max(X_min,int(ceil(X_frac*float(K[j]))))
                if K_lim[j] >= X:
                    mHG_n, mHG_s, pval = xlmhg.test(v,X,L,K=int(K[j]),mat=mat,pval_thresh=pval_thresh)

                    # check if GO term is significantly enriched
                    if pval <= pval_thresh:
                        # generate mHGTermResult
                        sel = np.nonzero(A[:,j])[0] # ranks of all the 1's
                        mHG_k_n = np.sum(sel < mHG_n) 
                        sel_genes = [ranked_genes[i] for i in sel]
                    
                        #def __init__(self,term,genes,pval,ranks,N,X,L,mHG_n=None,mHG_k_n=None,mHG_s=None):
                        #print terms[j],pval,sel_genes,sel
                        enr = mHGTermResult(terms[j],pval,sel_genes,sel,p,X,L,mHG_n,mHG_k_n,mHG_s)

                        enriched_terms.append(enr)

        self.message('done!')

        # calculate enrichment score
        self.message('Calculating enrichment score (using p-value threshold psi=%.1e) for enriched terms...', \
                escore_pval_thresh)
        for term in enriched_terms:
            term.calculate_escore(escore_pval_thresh)

        # report results
        q = len(enriched_terms)
        ignored = m - tested
        if ignored > 0:
            self.message('%d / %d GO terms (%.1f%%) had less than %d genes annotated with them and were ignored.', \
                        ignored,m,100*(ignored/float(m)),X_min)

        self.message('%d / %d tested GO terms were found to be significantly enriched (p-value <= %.1e).', \
                q,tested,pval_thresh)

        return enriched_terms
