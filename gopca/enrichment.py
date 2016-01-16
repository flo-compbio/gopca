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

"""Module for gene set enrichment (GSE) analysis using the XL-mHG test.

The `GSEAnalysis` class performs the tests, and the results are represented
by `GSEResult` objects.
"""

# TO-DO:
# * E-score calculation should be part of XL-mHG package
# * GSEResult should inherit from "mHGResult" (new class in XL-mHG package)

import sys
import os
import csv
import logging
import gzip
import cPickle as pickle
from math import ceil
from collections import Iterable

import numpy as np
from scipy.stats import hypergeom

from genometools import misc
from genometools.basic import GeneSet, GeneSetDB
from genometools.expression import ExpGenome

import xlmhg

logger = logging.getLogger(__name__)

class GSEResult(object):
    """Result of an XL-mHG-based test for gene set enrichment in a ranked list.

    Parameters
    ----------
    gene_set: `genometools.basic.GeneSet` object
        The gene set tested.
    stat: float
        The XL-mHG test statistic.
    pval: float
        The XL-mHG p-value.
    indices: list (or tuple / ndarray) of int
        The indices corresponding to the "1's" in the ranked list.
    genes: list or tuple of str
        The names of the genes corresponding to the "1's" in the ranked list.
   
    """

    #Note: Change this so that it inherits from class mHGResult?
    #     (only additional attributes: term, genes).

    # TO-DO: finish documentation

    def __init__(self, n, stat, pval, N, X, L,
            indices, gene_set, genes):

        assert isinstance(n, int)
        assert isinstance(stat, float)
        assert isinstance(pval, float)
        assert isinstance(N, int)
        assert isinstance(X, int)
        assert isinstance(L, int)
        assert isinstance(indices, Iterable)
        assert isinstance(gene_set, GeneSet)
        assert isinstance(genes, Iterable)

        self.n = n
        self.stat = stat
        self.pval = pval
        self.N = N
        self.X = X # XL-mHG "X" parameter
        self.L = L # XL-mHG "L" parameter

        self.indices = np.int32(indices).copy()
        self.indices.flags.writeable = False # makes it hashable

        self.gene_set = gene_set
        self.genes = tuple(genes)

        # strength of enrichment
        self.escore_pval_thresh = None
        self.escore = None

    def __repr__(self):
        return '<%s object (gene_set_id=%s; pval=%.1e; hash=%d)' \
                %(self.__class__.__name__, self.gene_set_id, hash(self))

    def __str__(self):
        return '<%s object (gene_set=%s; pval=%.1e)>' \
                %(self.__class__.__name__, str(self.gene_set), self.pval)

    def __hash__(self):
        data = []
        data.append(self.n)
        data.append(self.stat)
        data.append(self.pval)
        data.append(self.N)
        data.append(self.X)
        data.append(self.L)
        data.append(self.indices.data)
        data.append(self.gene_set)
        data.append(self.genes)
        return hash(tuple(data))

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) != type(other):
            return False
        else:
            return repr(self) == repr(other)

    def __ne__(self, other):
        return not (self == other)

    def __setstate__(self, d):
        self.__dict__ = d
        self.indices.flags.writeable = False

    @property
    def K(self):
        return self.indices.size

    @property
    def k_n(self):
        return int(np.sum(self.indices < self.n))

    def calculate_escore(self, pval_thresh):
        """ Calculate XL-mHG E-score.  """
        N = self.N
        K = self.K
        X = self.X
        L = self.L
        indices = self.indices
        
        if K == 0 or L == N or K < X:
            return 0

        k_max = 0
        fe_max = 0.0
        k = 1
        pval = 1.0

        while k <= K and indices[k-1] < L:
            if k >= X:
                n = indices[k-1] + 1
                if pval_thresh == 1.0 or hypergeom.sf(k-1,N,K,n) <= pval_thresh:
                    fe = k / (K * (n / float(N)))
                    if fe >= fe_max:
                        fe_max = fe
                        k_max = k
            k += 1

        self.escore_pval_thresh = pval_thresh
        self.escore = fe_max

    def get_pretty_format(self,omit_param=True,max_name_length=0):
        # TO-DO: clean up, commenting
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
        # TO-DO: clean up, commenting
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

class GSEAnalysis(object):
    """Test ranked gene lists for gene set enrichment using the XL-mHG test.

    The class is intialized with a set of valid gene names (an `ExpGenome`
    object), as well as a set of genesets (a `GeneSetDB` object). During
    initialization, a binary "gene-by-gene set" matrix is constructed,
    which stores information about which gene is contained in each gene set.
    This matrix is quite sparse, and requires a significant amount of memory.
    As an example, for a set of p = 10,000 genes and n = 10,000 gene sets, this
    matrix is of size 100 MB in the memory (i.e., p x n bytes).

    Once the class has been initialized, the function `get_enriched_gene_sets`
    can be called with a ranked list of genes, a significance threshold, and a
    set of parameters for the XL-mHG test. This function returns a list of
    `GSEResult` objects, one for each gene set that was found to be
    significantly enriched.

    Parameters
    ----------
    genome: `genometools.expression.ExpGenome` object
        See :attr:`genome` attribute.
    geneset_db: `genometools.basic.GeneSetDB` object
        See :attr:`geneset_db` attribute.

    Attributes
    ----------
    genome: `genometools.expression.ExpGenome` object
        The universe of genes.
    geneset_db: `genometools.basic.GeneSetDB` object
        The gene set database.
    """

    def __init__(self, genome, gene_set_db):

        assert isinstance(genome, ExpGenome)
        assert isinstance(gene_set_db, GeneSetDB)

        self.genome = genome
        self.gene_set_db = gene_set_db

        # generate annotation matrix by going over all gene sets
        logger.info('Generating gene x GO term matrix...')
        self.A = np.zeros((genome.p, gene_set_db.n), dtype = np.uint8)
        for j, gs in enumerate(self.gene_set_db.gene_sets):
            for g in gs.genes:
                try:
                    idx = self.genome.index(g)
                except ValueError:
                    pass
                else:
                    self.A[idx,j] = 1

    @property
    def genes(self):
        return self.genome.genes

    def get_enriched_gene_sets(self, ranked_genes, pval_thresh,
            X_frac, X_min, L,
            escore_pval_thresh = None, gene_set_ids = None, mat = None):
        """Tests gene set enrichment given a ranked list of genes.

        This function also calculates XL-mHG E-scores for the enriched gene
        sets, using ``escore_pval_thresh`` as the p-value threshold "psi".
        """

        # checks
        assert isinstance(ranked_genes, Iterable)
        for g in ranked_genes:
            assert isinstance(g, (str, unicode))
        assert isinstance(pval_thresh, float)
        assert isinstance(X_frac, float)
        assert isinstance(X_min, int)
        assert isinstance(L, int)

        if escore_pval_thresh is not None:
            assert isinstance(escore_pval_thresh, float)
        if gene_set_ids is not None:
            assert isinstance(gene_set_ids, Iterable)
            for id_ in gene_set_ids:
                assert isinstance(id_, (str, unicode))
        if mat is not None:
            assert isinstance(mat, np.ndarray)
        
        genes = self.genes
        gene_set_db = self.gene_set_db
        A = self.A

        if escore_pval_thresh is None:
            # if no separate E-score p-value threshold is specified, use the
            # p-value threshold (this results in very conservative E-scores)
            logger.warning('Setting the E-score p-value threshold to the ' +
            'global significance threshold results in conservative E-scores.')
            escore_pval_thresh = pval_thresh

        # test only some terms?
        if gene_set_ids is not None:
            gs_indices = np.int64([self.gene_set_db.index(id_)
                    for id_ in gene_set_ids])
            gene_sets = [gene_set_db[id_] for id_ in gene_set_ids]
            gene_set_db = GeneSetDB(gene_sets)
            A = A[:,gs_indices] # not a view!

        # reorder rows in annotation matrix to match the given gene ranking
        # also exclude genes not in the ranking
        gene_indices = np.int64([self.genome.index(g) for g in ranked_genes])
        A = A[gene_indices,:] # not a view either!

        # determine largest K
        K_lim = np.sum(A[:L,:], axis = 0, dtype = np.int64)
        K_rem = np.sum(A[L:,:], axis = 0, dtype = np.int64)
        K = K_lim + K_rem
        K_max = np.amax(K)

        # prepare matrix for XL-mHG p-value calculation
        p, m = A.shape
        if mat is None:
            mat = np.empty((K_max + 1, p + 1), dtype = np.longdouble)

        # find enriched GO terms
        logger.info('Testing %d gene sets for enrichment...', m)
        logger.debug('(N = %d, X_frac = %.2f, X_min = %d, L = %d; K_max = %d)', \
                    len(ranked_genes), X_frac, X_min, L, K_max)

        enriched = []
        tested = 0 # number of tests conducted
        N,m = A.shape
        for j in range(m):
            # determine gene set-specific value for X (based on K[j])
            X = max(X_min, int(ceil(X_frac * float(K[j]))))

            # determine significance of gene set enrichment using XL-mHG test
            # (only if there are at least X gene set genes in the list)
            if K[j] >= X:
                tested += 1

                # we only need to perform the XL-mHG test if there are enough
                # gene set genes at or above L'th rank (otherwise, pval = 1.0)
                if K_lim[j] >= X:
                    v = np.ascontiguousarray(A[:,j]) # copy
                    n, stat, pval = xlmhg.test(v, X, L, K = int(K[j]),
                            mat = mat, pval_thresh = pval_thresh)

                    # check if gene set is significantly enriched
                    if pval <= pval_thresh:
                        # generate GSEResult
                        sel = np.nonzero(A[:,j])[0] # indices of all the 1's
                        k_n = np.sum(sel < n) 
                        sel_genes = [ranked_genes[i] for i in sel]
                        result = GSEResult(n, stat, pval, N, X, L, sel, \
                                gene_set_db[j], sel_genes)
                        enriched.append(result)

        # calculate enrichment score
        logger.debug('Calculating enrichment score (using p-value threshold ' +
                'psi=%.1e) for enriched gene sets...', escore_pval_thresh)
        for result in enriched:
            result.calculate_escore(escore_pval_thresh)

        # report results
        q = len(enriched)
        ignored = m - tested
        if ignored > 0:
            logger.debug('%d / %d gene sets (%.1f%%) had less than X genes ' +
                    'annotated with them and were ignored.',
                    ignored,m,100*(ignored/float(m)))

        logger.info('%d / %d gene sets were found to be significantly ' +
                'enriched (p-value <= %.1e).', q, m, pval_thresh)

        return enriched
