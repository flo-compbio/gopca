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

import numpy as np

import gopca
from gopca.tools import misc
from gopca.xlmhg import xlmHG_cython as mHG
from gopca import enrichment_score

class mHGTermResult(object):
	"""
	Stores mHG result for one particular term.
	"""
	def __init__(self,term,p_value,N,n,K,genes,X=0,L=0):
		self.term = term
		self.p_value = p_value
		self.N = N
		self.n = n
		self.K = K
		self.genes = frozenset(genes)
		self.k = len(self.genes)
		self.X = X
		self.L = L
		if n == 0 or K == 0:
			self.fold_enrichment = float('nan')
		else:
			self.fold_enrichment = self.k/(n*(K/float(N)))

	def __repr__(self):
		return "<mHGTermResult: %s (p=%.1e; fe=%.1x; X=%d; L=%d; %d/%d@%d/%d), gene set %d>" \
				%(self.term.id,self.p_value,self.fold_enrichment,self.X,self.L,self.k,self.K,self.n,self.N,hash(self.genes))

	def __str__(self):
		return "<mHGTermResult of GO term '%s': p-value = %.1e, fold enrichment = %.2fx, %d/%d genes @ %d (N=%d)>" \
				%(str(self.term),self.p_value,self.fold_enrichment, self.k, self.K, self.n, self.N)

	def __hash__(self):
		return hash(repr(self))

	def __eq__(self,other):
		if type(self) != type(other):
			return False
		elif self.term == other.term and self.genes == other.genes:
			return True
		else:
			return False

	def get_pretty_format(self,GO=None,omit_acc=False,omit_param=True,nitty_gritty=True,max_name_length=0):
		term_str = '/'.join(self.term)
		if GO is not None:
			term = GO.terms[self.term[0]]
			term_str = term.get_pretty_format(omit_acc=omit_acc,max_name_length=max_name_length)
		details = ''
		if nitty_gritty:
			param_str = ''
			if not omit_param:
				param_str = ' (X=%d,L=%d)' %(self.X,self.L)
			details = ' [p=%.1e,e=%.2fx,%d/%d@%d%s]' %(self.p_value,self.fold_enrichment,len(self.genes),self.K,self.n,param_str)
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

	def test_enrichment(self,ranked_genes,pvalue_threshold,X,L,fold_enrichment_threshold=None,selected_terms=[],mat=None,quiet=False):

		genes = self.genes
		terms = self.terms
		term_ids = self.term_ids
		original_dimensions = self.A.shape
		A = self.A

		# test only some terms?
		if selected_terms:
			term_indices = np.int64([misc.bisect_index(term_ids,t) for t in selected_terms])
			A = A[:,term_indices] # not a view!

		# sort rows in annotation matrix (and exclude genes not in the ranking)
		order = []
		for g in ranked_genes:
			idx = misc.bisect_index(genes,g)
			order.append(idx)
		order = np.int64(order)
		A = A[order,:] # not a view either!

		# determine largest K
		K = np.sum(A,axis=0,dtype=np.int64)
		K_max = np.amax(K)

		# prepare matrix for XL-mHG p-value calculation
		p,m = A.shape
		if mat is None:
			mat = np.empty((K_max+1,p+1),dtype=np.longdouble)

		# find enriched GO terms
		if not quiet:
			print "Testing for enrichment of %d terms..." %(m)
			print "(N = %d, X = %d, L = %d; K_max = %d)" %(len(ranked_genes),X,L,K_max)
			sys.stdout.flush()

		enriched_terms = []
		tested = 0
		for j in range(m):
			#if j >= 1000: break
			if (not quiet) and (j % 100) == 0:
				print "\r%d..." %(j), ; sys.stdout.flush()

			v = np.ascontiguousarray(A[:,j]) # copy

			# determine significance of enrichment using XL-mHG test
			# (only if there are at least X genes annotated with this term)
			if K[j] >= X:
				tested += 1
				threshold,_,pval = mHG.mHG_test(v,p,K[j],L,X,mat,pvalue_threshold=pvalue_threshold)

				# check if GO term is significantly enriched
				if pval <= pvalue_threshold:
					sel = np.nonzero(A[:threshold,j])[0]
					sel_genes = [ranked_genes[i] for i in sel]
					k = sel.size

					enr = mHGTermResult(terms[j],pval,p,threshold,K[j],sel_genes)
					# test if associated fold enrichment exceeds fold enrichment threshold (if one was specified)
					if fold_enrichment_threshold is None or enr.fold_enrichment >= fold_enrichment_threshold:
						enriched_terms.append(enr)

		if not quiet:
			print 'done!'; sys.stdout.flush()
			if X > 0:
				ignored = m - tested
				print '%d/%d GO terms (%.1f%%) had less than %d genes annotated with them and were ignored.' %(ignored,m,100*(ignored/float(m)),X)
			q = len(enriched_terms)
			fold_term = ''
			if fold_enrichment_threshold is not None:
				fold_term = ' and fold enrichment >= %.1fx ' %(fold_enrichment_threshold)
			print '%d / %d tested GO terms (%.1f%%) were found to be enriched (p-value <= %.1e%s).' \
					%(q,tested,100*(q/float(tested)),pvalue_threshold,fold_term)
			sys.stdout.flush()

		assert self.A.shape == original_dimensions
		return enriched_terms

	def apply_thresholds(self,enrichments,pvalue_threshold,fold_enrichment_threshold=None,quiet=False):
		filtered = []

		terms = self.terms
		m = len(enrichments)
		for j,enr in enumerate(enrichments):
			if enr.p_value <= pvalue_threshold:
				if fold_enrichment_threshold is not None and enr.fold_enrichment < fold_enrichment_threshold:
					continue
				filtered.append(enr)

		if not quiet:
			kept = len(filtered)
		return filtered

	def get_enrichment_scores(self,enrichments,ranked_genes,X,L,p_max):

		genes = self.genes
		term_ids = self.term_ids
		A = self.A

		# select annotations matrix columns corresponding to enriched GO terms
		enriched_terms = [enr.term[0] for enr in enrichments]
		term_indices = np.int64([misc.bisect_index(term_ids,t) for t in enriched_terms])
		A = A[:,term_indices] # not a view!

		# sort rows in annotation matrix (and exclude genes not in the ranking)
		indices = []
		for g in ranked_genes:
			idx = misc.bisect_index(genes,g)
			indices.append(idx)
		indices = np.int64(indices)
		A = A[indices,:] # not a view either!

		q = len(enrichments)
		es = np.zeros(q,dtype=np.float64)
		p = len(ranked_genes)
		K = np.sum(A,axis=0,dtype=np.int64)
		for i in range(q):
			es[i] = enrichment_score.get_enrichment_score(A[:,i],p,K[i],X,L,p_max)
		return es

	def pretty_print_enrichment(self):
		raise NotImplemented
		N,M = self.A.shape
		pvalues = np.empty(M,dtype=np.float64)
