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
from sklearn.decomposition import PCA

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
	parser.add_argument('-X','--go-mHG-X',type=int,default=5) # 0=off
	parser.add_argument('-L','--go-mHG-L',type=int,default=1000) #0=off

	### options for selecting the number of PCs to test

	# direct method:
	# --num-pc directly specifies the number of PCs to test
	# --pc-max forces testing to stop after the first X PCs (for use in combination with other methods)
	parser.add_argument('--pc-num',type=int,default=0) # 0=off
	parser.add_argument('--pc-max',type=int,default=15) # 0=off

	# data-driven method:
	# stop testing PCs when the cumulative variance explained surpasses X %
	parser.add_argument('--pc-cum-var',type=float,default=80.0) # in percent (0=off)

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

	# exclude genes from most enriched term
	new_ranked_genes = []
	new_L = L
	for i,g in enumerate(ranked_genes):
		if g not in genes_used:
			new_ranked_genes.append(g)
		elif i < L: # gene was already used, adjust L if necessary
			new_L -= 1
	ranked_genes = new_ranked_genes
	L = new_L

	# filter enrichments
	K_max = max([enr.K for enr in todo])
	p = len(ranked_genes)
	mat = np.zeros((K_max+1,p+1),dtype=np.longdouble)
	while todo:
		most_enriched = todo[0]
		term_id = most_enriched.term[0]

		# test if enrichment is still significant after removing all previously used genes
		enr = M.test_enrichment(new_ranked_genes,pvalue_threshold,X,new_L,selected_terms=[term_id],mat=mat,quiet=True)[0]
		if enr.p_value <= pvalue_threshold and enr.fold_enrichment >= fold_enrichment_threshold:
			# if so, keep it!
			kept.append(most_enriched)
			# next, exclude selected genes from further analysis: 1) adjust L 2) update set of excluded genes
			genes_used.update(most_enriched.genes) # add selected genes to set of used genes
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
	annotation_files = [args.annotation_gene_file,args.annotation_term_file,args.annotation_matrix_file]

	# parameters
	go_pvalue_threshold = args.go_pvalue_threshold
	go_fold_enrichment_threshold = args.go_fold_enrichment_threshold
	go_mHG_X = args.go_mHG_X
	go_mHG_L = args.go_mHG_L

	pc_num = args.pc_num
	pc_max = args.pc_max
	pc_cum_var = args.pc_cum_var/100.0

	# parameter checks?

	# misc options
	verbose = args.verbose

	# read expression data
	genes,samples,E = common.read_expression(expression_file)
	print "Expression matrix dimensions:", E.shape; sys.stdout.flush()

	if go_mHG_L == 0:
		go_mHG_L = len(genes)

	# read GO data
	print "Loading GO annotations...", ; sys.stdout.flush()
	GO = pickle.load(open(args.go_pickle_file))
	print "done!"; sys.stdout.flush()

	# create GOEnrichment object
	print "Reading GO annotation data...", ; sys.stdout.flush()
	M_enrich = GOEnrichment()
	M_enrich.read_annotations(*annotation_files)
	m = len(M_enrich.terms)
	print "read data for %d GO terms!" %(m); sys.stdout.flush()

	# determine number of PCs to compute
	compute_pc = pc_num
	if pc_num == 0:
		if pc_max > 0:
			compute_pc = pc_max
		else:
			compute_pc = 15 # if nothing is specified, look at first 15 PCs

	# there are at most n principal components (n = # samples)
	compute_pc = min(compute_pc,len(samples))
	# there are at most p-1 principal components (p = # genes)
	compute_pc = min(compute_pc,len(genes)-1)

	# perform PCA
	print 'Performing PCA...', ; sys.stdout.flush()
	sys.stdout.flush()
	M_pca = PCA(n_components = compute_pc)
	M_pca.fit(E.T)
	print 'done!'
	frac = M_pca.explained_variance_ratio_
	cum_frac = np.cumsum(frac)
	print "Cumulative fraction of variance explained for first %d PCs:" %(compute_pc)
	print cum_frac
	sys.stdout.flush()

	# determine the number of PCs to test
	test_pc = compute_pc
	if pc_num == 0:
		# apply variance criterion
		if pc_cum_var > 0:
			sel = np.nonzero(cum_frac > pc_cum_var)[0]
			if sel.size > 0:
				test_pc = sel[0] + 1

	assert test_pc <= compute_pc
	print "GO-PCA will test the first %d principal components!" %(test_pc)
	sys.stdout.flush()

	# run GO-PCA!
	W = M_pca.components_.T
	print W.shape
	final_enrichment = []
	p = len(genes)
	res_var = None
	all_genes = set(genes)
	total_var = 0.0
	for pc in range(test_pc):

		print
		print '-------------------------------------------------------------------------'
		print "PC %d explains %.1f%% of the total variance." %(pc+1,100*frac[pc])
		total_var += frac[pc]
		print "The new cumulative fraction of total variance explained is %.1f%%." %(100*total_var)
		sys.stdout.flush()

		print "Testing for GO enrichment...", ; sys.stdout.flush()
		enrichment_dsc = get_pc_go_enrichment(M_enrich,pc+1,-W[:,pc],genes,go_mHG_X,go_mHG_L,go_pvalue_threshold,go_fold_enrichment_threshold)
		enrichment_asc = get_pc_go_enrichment(M_enrich,-pc-1,W[:,pc],genes,go_mHG_X,go_mHG_L,go_pvalue_threshold,go_fold_enrichment_threshold)
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
