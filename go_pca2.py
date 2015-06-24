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

# allow explicit relative imports in executable script
# source: http://stackoverflow.com/a/6655098
if __name__ == '__main__' and __package__ is None:
	parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	sys.path.insert(1, parent_dir)
	import gopca
	__package__ = 'gopca'

import argparse
import csv
import cPickle as pickle
import itertools as it

import numpy as np
from scipy import stats
from sklearn.decomposition import RandomizedPCA

from gopca.tools import misc
from gopca import common
from gopca.go_enrichment import GOEnrichment
from gopca import fdr

#from goparser.parser import GOParser
from .go_pca_objects import mHGTermResultWithPC,SignatureMatrix,NamedGeneSet

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	###
	### Required arguments
	###

	# input files
	parser.add_argument('-e','--expression-file',required=True)
	parser.add_argument('-g','--go-pickle-file',required=True)

	parser.add_argument('-ag','--annotation-gene-file',required=True)
	parser.add_argument('-at','--annotation-term-file',required=True)
	parser.add_argument('-am','--annotation-matrix-file',required=True)

	# output file
	parser.add_argument('-o','--output-file',required=True)

	###
	### Optional arguments
	###

	# main parameters
	parser.add_argument('-p','--go-pvalue-threshold',type=float,default=1e-6)
	parser.add_argument('-f','--go-fold-enrichment-threshold',type=float,default=1.5)
	parser.add_argument('-s','--seed',type=int,default=-1)

	# options for testing associations between genes and PCs
	parser.add_argument('-z','--gene-pc-zscore-file',default=None)
	parser.add_argument('--num-pc',type=int,default=0) # pre-determined (0=off)

	parser.add_argument('--num-pc-default',type=int,default=15) # only if z-score file is provided
	parser.add_argument('--pc-fdr',type=float,default=0.01) # only if z-score file is provided
	parser.add_argument('--pc-min-genes',type=int,default=100) # only if z-score file is provided (0=off)
	parser.add_argument('--pc-max',type=int,default=0) # hard limit (0=off)
	parser.add_argument('--pc-min-variance',type=float,default=1.0,help='in percent') # ignore PCs capturing less than this fraction of variance (0 = off)

	parser.add_argument('--mHG-X',type=int,default=5)
	parser.add_argument('--mHG-L',type=int,default=1000)
	parser.add_argument('--enrichment-fdr',type=float,default=0.05)

	# output verbosity
	#parser.add_argument('-v','--verbose',action='store_true')
	parser.add_argument('-v','--verbose',action='store_false')

	return parser.parse_args()

def print_enrichment(l,GO):
	a = sorted(range(len(l)),key=lambda i: -l[i].fold_enrichment)
	for i in a:
		enr = l[i]
		term = GO.terms[enr.term[0]]
		goterm_genes = GO.get_goterm_genes(term.id)
		print enr.get_pretty_format(GO)

def get_pc_go_enrichment(M,pc,w,genes,X,L,pvalue_threshold,fold_enrichment_threshold,quiet=False):
	a = np.argsort(w)
	ranked_genes = [genes[i] for i in a]
	# test enrichment
	enrichment = get_go_enrichment(M,ranked_genes,X,L,pvalue_threshold,fold_enrichment_threshold,quiet=quiet)
	# add PC information to enrichment results
	enrichment = [mHGTermResultWithPC.from_mHGTermResult(pc,enr) for enr in enrichment]
	return enrichment

def get_go_enrichment(M,ranked_genes,X,L,pvalue_threshold,fold_enrichment_threshold,filtering=True,quiet=False):
	# tests for enrichment, applies thresholds, performs filtering to reduce redundancy
	enrichment = M.test_enrichment(ranked_genes,pvalue_threshold,X,L,quiet=quiet)
	enrichment = M.apply_thresholds(enrichment,pvalue_threshold,fold_enrichment_threshold,quiet=quiet)
	if filtering:
		enrichment = filter_enriched_go_terms(enrichment,M,ranked_genes,X,L,pvalue_threshold,fold_enrichment_threshold,quiet=quiet)
	return enrichment

def filter_enriched_go_terms(enrichment,M,ranked_genes,X,L,pvalue_threshold,fold_enrichment_threshold,quiet=False):

	if len(enrichment) <= 1:
		return []

	# sort enrichments by fold change
	todo = sorted(enrichment,key=lambda enr: -enr.fold_enrichment)
	most_enriched = todo[0]
	genes_used = set(most_enriched.genes)
	kept = [most_enriched]
	todo = todo[1:]
	ranked_genes = ranked_genes[:] # make a copy here!

	# filter enrichments
	K_max = max([enr.K for enr in todo])
	p = len(ranked_genes)
	mat = np.zeros((K_max+1,p+1),dtype=np.longdouble)
	while todo:
		most_enriched = todo[0]
		term_id = most_enriched.term[0]
		ranked_genes = [g for g in ranked_genes if g not in genes_used]

		# test if enrichment is still significant after removing all previously used genes
		enr = M.test_enrichment(ranked_genes,pvalue_threshold,X,L,selected_terms=[term_id],mat=mat,quiet=True)[0]
		if enr.p_value <= pvalue_threshold and enr.fold_enrichment >= fold_enrichment_threshold:
			# if so, keep it!
			kept.append(most_enriched)
			# next, exclude selected genes from further analysis: 1) adjust L 2) update set of excluded genes
			# 1) figure out how many genes above L were selected, and adjust L accordingly
			exclude_genes = genes_used - set(most_enriched.genes) # newly excluded genes
			new_L = L
			for i in range(L):
				if ranked_genes[i] in exclude_genes:
					new_L -= 1
			L = new_L
			genes_used.update(most_enriched.genes) # add selected genes to set of used genes

		todo = todo[1:]

	if not quiet:
		print 'Filtering: kept %d / %d enrichments.' %(len(kept),len(enrichment)); sys.stdout.flush()

	return kept

def remove_redundant_terms(new_enrichment,previous_enrichment,GO):
	if len(previous_enrichment) == 0:
		return new_enrichment
	keep = []
	previous_terms = set([enr.term[0] for enr in previous_enrichment])
	for enr in new_enrichment:
		term_id = enr.term[0]
		term = GO.terms[term_id]
		novel = True
		for t in set([term_id]) | term.ancestors | term.descendants:
			if t in previous_terms:
				novel = False
				break
		if novel:
			keep.append(enr)
	return keep

def main(args):

	# input files
	expression_file = args.expression_file
	zscore_file = args.gene_pc_zscore_file
	annotation_files = [args.annotation_gene_file,args.annotation_term_file,args.annotation_matrix_file]

	# parameters
	go_pvalue_threshold = args.go_pvalue_threshold
	go_fold_enrichment_threshold = args.go_fold_enrichment_threshold
	pc_fdr = args.pc_fdr
	pc_min_genes = args.pc_min_genes
	pc_min_var = args.pc_min_variance/100.0
	pc_max = args.pc_max
	enrichment_fdr = args.enrichment_fdr
	mHG_X = args.mHG_X
	mHG_L = args.mHG_L
	num_pc = args.num_pc
	seed = args.seed

	# parameter checks?

	# select/generate seed for random number generator
	max_int = np.iinfo(np.int32).max
	if seed < 0:
		seed = np.random.randint(0,max_int)
	np.random.seed(seed)
	print "Seed used:",seed; sys.stdout.flush()

	#workers = args.workers

	# misc options
	verbose = args.verbose

	# read expression data
	genes,samples,E = common.read_expression(expression_file)
	print "Expression matrix dimensions:",E.shape

	# read GO data
	print "Loading GO annotations...", ; sys.stdout.flush()
	GO = pickle.load(open(args.go_pickle_file))
	print "done!"; sys.stdout.flush()

	# get number of GO terms to be tested
	m = None
	with open(args.annotation_term_file) as fh:
		m = len(fh.readlines())
	print "Testing %d GO terms." %(m)

	# create GOEnrichment object
	print "Preparing GO enrichment tests...", ; sys.stdout.flush()
	M = GOEnrichment()
	M.read_annotations(*annotation_files)
	print "done!"; sys.stdout.flush()

	# load pc/gene z-scores
	S = None
	Z = None
	k = None
	crit_p = None
	if zscore_file is not None:
		_,_,Z = common.read_gene_data(zscore_file)
		#print E.shape[0],Z.shape[0]
		assert Z.shape[0] == E.shape[0] # make sure number of genes is correct

		# convert Z-scores to P-values (= one-sided tail)
		S = stats.norm.sf(Z,loc=0.0,scale=1.0)
		# calculate FDR thresholds
		fdr_results = [fdr.fdr_bh(pval,pc_fdr) for pval in S.T]
		k = np.int64([f[0] for f in fdr_results])
		print "Significant genes per PC:"
		print k
		sys.stdout.flush()

	# determine number of PCs to compute
	if num_pc == 0:
		if S is not None:
			if pc_min_genes > 0:
				sel = np.nonzero(k < pc_min_genes)[0]
				if sel.size > 0:
					num_pc = sel[0]
				else:
					num_pc = S.shape[1]
			else:
				num_pc = S.shape[1]
		else:
			num_pc = num_pc_default
		if pc_max > 0:
			num_pc = min(num_pc,pc_max)

	print "Will compute the first %d principal components." %(num_pc)
	sys.stdout.flush()

	# calculate PCA and determine number of PCs to test
	P = RandomizedPCA(n_components = num_pc)
	P.fit(E.T)
	frac = P.explained_variance_ratio_
	print "Fraction variance explained per PC:"
	print frac
	print np.cumsum(frac)
	sys.stdout.flush()

	# apply variance criterion
	if pc_min_var > 0:
		sel = np.nonzero(frac < pc_min_var)[0]
		if sel.size > 0:
			num_pc = min(num_pc,sel[0])

	print "Will test the first d=%d principal components!" %(num_pc)
	sys.stdout.flush()

	# run GO-PCA!
	C = P.components_.T
	final_enrichment = []
	p = len(genes)
	res_var = None
	all_genes = set(genes)
	total_var = 0.0
	d = C.shape[1]
	for pc in range(num_pc):

		sig_pos = mHG_L
		sig_neg = mHG_L
		if S is not None and enrichment_fdr > 0:
			# we have Z-scores, and user wants us to use them to help set L for each PC
			_,crit_p = fdr.fdr_bh(S[:,pc],enrichment_fdr)
			sel = np.zeros(p,dtype=np.bool_)
			sel[np.nonzero(S[:,pc] <= crit_p)[0]] = True
			sig_pos = min(sig_pos,np.sum(np.all(np.c_[sel,Z[:,pc]>0],axis=1)))
			sig_neg = min(sig_neg,np.sum(np.all(np.c_[sel,Z[:,pc]<0],axis=1)))
			
		s = ''
		if S is not None:
			s = 'has %d significant genes (at FDR=%.2f) and ' %(k[pc],pc_fdr)
		print
		print '-------------------------------------------------------------------------'
		print "PC %d %sexplains %.1f%% of the total variance." %(pc+1,s,100*frac[pc])
		total_var += frac[pc]
		print "The new cumulative fraction of total variance explained is %.1f%%." %(100*total_var)
		sys.stdout.flush()

		print "Testing for GO enrichment...", ; sys.stdout.flush()
		enrichment_dsc = get_pc_go_enrichment(M,pc+1,-C[:,pc],genes,mHG_X,mHG_L,go_pvalue_threshold,go_fold_enrichment_threshold)
		enrichment_asc = get_pc_go_enrichment(M,-pc-1,C[:,pc],genes,mHG_X,mHG_L,go_pvalue_threshold,go_fold_enrichment_threshold)
		enrichment = enrichment_dsc + enrichment_asc

		print "# enriched terms:",len(enrichment); sys.stdout.flush()
		before = len(enrichment)

		enrichment = remove_redundant_terms(enrichment,final_enrichment,GO)
		removed = before - len(enrichment)
		print "Kept %d non-redundant enrichments (%d removed)." %(len(enrichment),removed)
	
		print_enrichment(enrichment,GO)
		final_enrichment.extend(enrichment)
		print "# total enriched terms:", len(final_enrichment); sys.stdout.flush()

		pc += 1

	print_enrichment(final_enrichment,GO)

	with open(args.output_file,'w') as ofh:
		pickle.dump(final_enrichment,ofh,pickle.HIGHEST_PROTOCOL)

	print "Done!"
	return 0

if __name__ == '__main__':
	return_code = main(read_args_from_cmdline())
	sys.exit(return_code)
