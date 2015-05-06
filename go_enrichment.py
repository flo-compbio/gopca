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

from tools import misc
from mhg import xlmHG_cython as mHG

class mHGTermResult(object):
	"""
	Stores mHG result for one particular term.
	"""
	def __init__(self,term,p_value,N,n,K,genes):
		self.term = term
		self.p_value = p_value
		self.N = N
		self.n = n
		self.K = K
		genes = frozenset(genes)
		self.k = len(genes)
		self.genes = genes
		if n == 0 or K == 0:
			self.fold_enrichment = float('nan')
		else:
			self.fold_enrichment = self.k/(n*(K/float(N)))

	def __repr__(self):
		return "<mHGTermEnrichment of term '%s', %d genes (hash:)>" %(self.term.id,self.k,hash(self.genes))

	def __str__(self):
		return "<mHG_Enrichment of term '%s': p-value = %.1e, fold enrichment = %.2fx, %d/%d genes @ %d>" \
				%(str(self.term),self.p_value,self.fold_enrichment, self.k, self.K, self.n)

	def __hash__(self):
		return hash(repr(self))

	def __eq__(self,other):
		if type(self) != type(other):
			return False
		elif self.term == other.term and self.genes == other.genes:
			return True
		else:
			return False

	def get_pretty_format(self,GO=None,omit_acc=False,nitty_gritty=True,max_name_length=0):
		term_str = '/'.join(self.term)
		if GO is not None:
			term = GO.terms[self.term[0]]
			term_str = term.get_pretty_format(omit_acc=omit_acc,max_name_length=max_name_length)
		details = ''
		if nitty_gritty:
			details = ' [p=%.1e,e=%.2fx,%d/%d@%d]' %(self.p_value,self.fold_enrichment,len(self.genes),self.K,self.n)
		return '%s%s' %(term_str,details)


class GOEnrichment(object):

	def __init__(self):
		self.genes = None
		self.terms = None
		self.A = None

	def read_annotations(self,annotation_gene_file,annotation_term_file,annotation_matrix_file):
		# read annotation data
		self.genes = misc.read_single(annotation_gene_file) # we assume genes are sorted alphabetically
		self.terms = misc.read_all(annotation_term_file)
		self.term_ids = [t[0] for t in self.terms]
		self.A = np.load(gzip.open(annotation_matrix_file))
		N,M = self.A.shape
		assert len(self.genes) == N
		assert len(self.terms) == M

	def get_term_index(self,term_id):
		return misc.bisect_index(self.term_ids,term_id)

	def test_enrichment(self,ranked_genes,X,L,selected_terms=None,quiet=False):

		enrichments = []

		genes = self.genes
		terms = self.terms
		A = self.A.copy()

		# test only some terms?
		if selected_terms is not None:
			assert len(selected_terms) > 0
			selected_terms = np.int64(selected_terms)
			terms = [terms[j] for j in selected_terms]
			A = A[:,selected_terms]

		# sort rows in annotation matrix (and exclude genes not in the ranking)
		order = []
		for g in ranked_genes:
			idx = misc.bisect_index(genes,g)
			order.append(idx)
		order = np.int64(order)
		A = A[order,:]

		N,M = A.shape
		K = np.sum(A,axis=0,dtype=np.int64)
		K_max = np.amax(K)
		if not quiet: print "K_max:",K_max
		

		pval = np.ones(M,dtype=np.float64)
		n = np.zeros(M,dtype=np.int64)
		k = np.zeros(M,dtype=np.int64)
		mat = np.empty((K_max+1,N+1),dtype=np.longdouble)

		if not quiet: print "Testing for enrichment of %d terms..." %(M), ; sys.stdout.flush()
		tests = 0
		v = np.empty(N,dtype=np.uint8)
		for j in range(M):
			#if j >= 1000: break
			if (j+1) % 100 == 0:
				if not quiet: print "%d..." %(j+1), ; sys.stdout.flush()

			v[:] = A[:,j] # copy

			# mHG
			if K[j] >= X:
				threshold,_,p = mHG.mHG_test(v,N,int(K[j]),L,X,mat)
			else:
				threshold = 0
				p = 1.0
			pval[j] = p
			n[j] = threshold
			if threshold > 0:
				k[j] = np.sum(v[:threshold],dtype=np.int64)
			sel = np.nonzero(A[:threshold,j])[0]
			sel_genes = [ranked_genes[i] for i in sel]
			assert k[j] == len(sel_genes)
			enr = mHGTermResult(terms[j],p,N,threshold,K[j],sel_genes)
			enrichments.append(enr)

		if not quiet: print "done!"; sys.stdout.flush()
		return enrichments

	def apply_thresholds(self,enrichments,p_value_threshold,fold_enrichment_threshold=None,quiet=False):
		filtered = []

		terms = self.terms
		m = len(enrichments)
		for j,enr in enumerate(enrichments):
			if enr.p_value <= p_value_threshold:
				if fold_enrichment_threshold is not None and enr.fold_enrichment < fold_enrichment_threshold:
					continue
				filtered.append(enr)

		if not quiet:
			kept = len(filtered)
			fold_term = ''
			if fold_enrichment_threshold is not None:
				fold_term = 'and fold enrichment >= %.1f ' %(fold_enrichment_threshold)
			print '%d / %d terms with p-value <= %.1e %s(%.1f%%)' %(kept,m,p_value_threshold,fold_term,100*(kept/float(m)))
		return filtered

	def prune_enrichment(self,GO):

		assert self.enrichment_filtered is not None
		self.enrichment_pruned = {}

		genes = self.genes
		significant_terms = sorted(self.enrichment_filtered.keys())
		m = len(significant_terms)

		filtered = np.zeros(m,dtype=np.bool_)
		for j,term_id in enumerate(significant_terms):
			enr = self.enrichment_filtered[term_id]
			fe = enr.fold_enrichment
			term = GO.terms[term_id]

			# filter rules: remove if
			# 	1) any descendent term (is also significant and) has greater or equal fold enrichment
			#	2) any ancestor term (is also significant and) has greater fold enrichment
			idx = None
			remove = False
			# go over ancestors
			for t in term.ancestors:
				try:
					if self.enrichment_filtered[t].fold_enrichment > fe:
						remove = True
						break
				except KeyError:
					pass
			if remove:
				filtered[j] = True
				break
		
			# go over descendants
			for t in term.descendants:
				try:
					if self.enrichment_filtered[t].fold_enrichment >= fe:
						remove = True
						break
				except KeyError:
					pass
			if remove:
				filtered[j] = True

		keep = np.nonzero(np.invert(filtered))[0]
		filtered_terms = [significant_terms[j] for j in keep]
		for t in filtered_terms:
			self.enrichment_pruned[t] = self.enrichment_filtered[t]
		print '%d / %d terms after pruning (%.1f%%)' %(keep.size,m,100*(keep.size/float(m))); sys.stdout.flush()

	def pretty_print_enrichment(self):
		raise NotImplemented
		N,M = self.A.shape
		pvalues = np.empty(M,dtype=np.float64)
