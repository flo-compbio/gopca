#!/usr/bin/env python

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
import gzip
import cPickle as pickle
from math import ceil

import numpy as np
from scipy.stats import hypergeom

import gopca
from gopca.tools import misc
from gopca.xlmhg.mHG import mHG_test
#from gopca import enrichment_score

class mHGTermResult(object):
	"""
	Stores mHG result for one particular term.

	Note: Change this so that it inherits from class mHGTerm
		  (only additional attributes: term, genes).
	"""

	def __init__(self,term,genes,pval,ranks,N,X,L,mHG_n=None,mHG_k_n=None,mHG_s=None):

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

	def __repr__(self):
		return "<mHGTermResult: %s (pval=%.1e; X=%d; L=%d; N=%d); genes hash=%d; positions hash=%d)>" \
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

	def get_max_fold_enrichment(self,pval_thresh=1.0):
		""" Find maximum fold enrichment. Returns both the number of genes selected and the corresponding fold enrichment. """
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

		return k_max, fe_max

	def get_pretty_format(self,omit_param=True,max_name_length=0):
		term_name = self.term[3]
		if max_name_length > 0 and len(term_str) > max_name_length:
			assert max_name_length >= 3
			term_name = term_str[:(len(term_str)-3)] + '...'
		term_str = term_name + ' (%d)' %(len(self.genes))
		param_str = ''
		if not omit_param:
			param_str = ' [X=%d,L=%d,N=%d]' %(self.X,self.L,self.N)
		details = ', p-value=%.1e%s' %(self.pval,param_str)
		return '%s%s' %(term_str,details)
		
	def get_pretty_GO_format(self,GO,omit_acc=False,omit_param=True,max_name_length=0):
		# accepts a GOParser object ("GO")
		term = GO.terms[self.term[0]]
		term_name = term.get_pretty_format(omit_acc=omit_acc,max_name_length=max_name_length)
		term_str = term_name + ' (%d)' %(len(self.genes))
		param_str = ''
		if not omit_param:
			param_str = ' [X=%d,L=%d,N=%d]' %(self.X,self.L,self.N)
		details = ', p-value=%.1e%s' %(self.pval,param_str)
		return '%s%s' %(term_str,details)

class GOEnrichment(object):

	def __init__(self):
		self.genes = None
		self.terms = None
		self.A = None

	def read_annotations(self,annotation_gene_file,annotation_term_file,annotation_matrix_file):
		# read annotation data
		self.genes = misc.read_single(annotation_gene_file) # we assume genes are sorted alphabetically
		self.terms = misc.read_all(annotation_term_file) # we assume terms are sorted alphabetically by term ID
		self.term_ids = [t[0] for t in self.terms]
		self.A = np.load(gzip.open(annotation_matrix_file))
		p,m = self.A.shape
		assert len(self.genes) == p
		assert len(self.terms) == m

	def test_enrichment(self,ranked_genes,pval_thresh,X_frac,X_min,L,selected_term_ids=[],mat=None,quiet=False):
		"""
		Tests GO term enrichment of either all terms or the terms specified by ``selected_term_ids''.
		"""

		genes = self.genes
		terms = self.terms
		term_ids = self.term_ids
		original_dimensions = self.A.shape
		A = self.A

		# test only some terms?
		if selected_term_ids:
			term_indices = np.int64([misc.bisect_index(term_ids,t) for t in selected_term_ids])
			terms = [self.terms[i] for i in term_indices]
			A = A[:,term_indices] # not a view!

		# sort rows in annotation matrix (and exclude genes not in the ranking)
		order = []
		for g in ranked_genes:
			idx = misc.bisect_index(genes,g)
			order.append(idx)
		order = np.int64(order)
		A = A[order,:] # not a view either!

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
		if not quiet:
			print "Testing for enrichment of %d terms..." %(m)
			print "(N = %d, X_frac = %.2f, X_min = %d, L = %d; K_max = %d)" \
					%(len(ranked_genes),X_frac,X_min,L,K_max)
			sys.stdout.flush()

		enriched_terms = []
		tested = 0
		tested_mHG = 0
		for j in range(m):
			#if j >= 1000: break
			if (not quiet) and (j % 100) == 0:
				print "\r%d..." %(j), ; sys.stdout.flush()

			v = np.ascontiguousarray(A[:,j]) # copy

			# determine significance of enrichment using XL-mHG test
			# (only if there are at least X_min genes annotated with this term before the L'th threshold)
			#if K_lim[j] >= X_min:
			if K[j] >= X_min:
				tested += 1
				# determine term-specific X (based on K[j])
				X = max(X_min,int(ceil(X_frac*float(K[j]))))
				if K_lim[j] >= X:
					tested_mHG += 1
					mHG_n, mHG_s, pval = mHG_test(v,X,L,K=int(K[j]),mat=mat,pval_thresh=pval_thresh)

					# check if GO term is significantly enriched
					if pval <= pval_thresh:
						# generate mHGTermResult
						sel = np.nonzero(A[:,j])[0] # ranks of all the 1's
						mHG_k_n = np.sum(sel < mHG_n) 
						sel_genes = [ranked_genes[i] for i in sel]
					
						#def __init__(self,term,genes,pval,ranks,N,X,L,mHG_n=None,mHG_k_n=None,mHG_s=None):
						enr = mHGTermResult(terms[j],sel_genes,pval,sel,p,X,L,mHG_n,mHG_k_n,mHG_s)

						enriched_terms.append(enr)

		if not quiet:
			print 'done!'; sys.stdout.flush()
			remaining = ''
			if X_min > 0:
				ignored = m - tested
				remaining = ' remaining'
				print ' %d/ %d GO terms (%.1f%%) had less than %d genes annotated with them and were ignored.' \
						%(ignored,m,100*(ignored/float(m)),X_min)
			if X_min > 0 or X_frac > 0:
				print 'Ran XL-mHG test on %d / %d%s GO terms (%.1f%%)' \
						%(tested_mHG,tested,remaining,100*(tested_mHG/float(tested)))
			q = len(enriched_terms)
			print '%d / %d tested GO terms (%.1f%%) were found to be significantly enriched (p-value <= %.1e).' \
					%(q,tested,100*(q/float(tested)),pval_thresh)
			sys.stdout.flush()

		assert self.A.shape == original_dimensions
		return enriched_terms

	def pretty_print_enrichment(self):
		raise NotImplemented
		N,M = self.A.shape
		pvalues = np.empty(M,dtype=np.float64)
