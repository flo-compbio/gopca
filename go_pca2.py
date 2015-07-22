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
from .go_pca_objects import GOPCAResult,GOPCASignature

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
	parser.add_argument('-L','--go-mHG-L',type=int,default=1000) # 0=off

	### options for selecting the number of PCs to test

	# direct method:
	# --pc-numc directly specifies the number of PCs to test
	# --pc-max forces testing to stop after the first X PCs (can be used in combination with --pc-cum-var)
	parser.add_argument('--pc-num',type=int,default=0) # 0=off
	parser.add_argument('--pc-max',type=int,default=15) # 0=off

	# data-driven method:
	# stop testing PCs when the cumulative variance explained surpasses X %
	parser.add_argument('--pc-cum-var',type=float,default=80.0) # in percent, 0=off

	# output verbosity
	#parser.add_argument('-v','--verbose',action='store_true')
	parser.add_argument('-v','--verbose',action='store_false')

	return parser.parse_args()

def print_signatures(signatures,GO):
	a = sorted(range(len(signatures)),key=lambda i: -signatures[i].fold_enrichment)
	for i in a:
		sig = signatures[i]
		term = GO.terms[sig.term[0]]
		goterm_genes = GO.get_goterm_genes(term.id)
		print sig.get_pretty_format(GO)

def get_pc_signatures(M,W,pc,genes,X,L,pvalue_threshold,fold_enrichment_threshold,quiet=False):
	pc_index = abs(pc)-1
	a = np.argsort(W[:,pc_index])
	if pc > 0:
		a = a[::-1]
	ranked_genes = [genes[i] for i in a]
	# get signatures
	signatures = get_signatures(M,ranked_genes,X,L,pvalue_threshold,fold_enrichment_threshold,quiet=quiet)
	# add PC information to signatures
	signatures = [GOPCASignature.from_mHGTermResult(pc,sig) for sig in signatures]
	return signatures

def get_signatures(M,ranked_genes,X,L,pvalue_threshold,fold_enrichment_threshold,filtering=True,quiet=False):
	# tests for enrichment, applies thresholds, performs filtering to reduce redundancy
	signatures = M.test_enrichment(ranked_genes,pvalue_threshold,X,L,quiet=quiet)
	signatures = M.apply_thresholds(signatures,pvalue_threshold,fold_enrichment_threshold,quiet=quiet)
	if filtering:
		signatures = filter_signatures(signatures,M,ranked_genes,X,L,pvalue_threshold,fold_enrichment_threshold,quiet=quiet)
	return signatures

def filter_signatures(signatures,M,ranked_genes,X,L,pvalue_threshold,fold_enrichment_threshold,quiet=False):

	if len(signatures) <= 1:
		return signatures

	# sort signatures by fold enrichment
	todo = sorted(signatures,key=lambda enr: -enr.fold_enrichment)
	most_enriched_signature = todo[0]
	genes_used = set(most_enriched_signature.genes)
	kept_signatures = [most_enriched_signature]
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

	# filter signatures
	K_max = max([sig.K for sig in todo])
	p = len(ranked_genes)
	mat = np.zeros((K_max+1,p+1),dtype=np.longdouble)
	while todo:
		#print '.', ; sys.stdout.flush()
		most_enriched_signature = todo[0]
		term_id = most_enriched_signature.term[0]

		# test if enrichment underlying signature is still significant after removing all previously used genes
		enr = M.test_enrichment(ranked_genes,pvalue_threshold,X,new_L,selected_terms=[term_id],mat=mat,quiet=True)
		enr = enr[0]
		#print '/', ; sys.stdout.flush()
		if enr.p_value <= pvalue_threshold and enr.fold_enrichment >= fold_enrichment_threshold:
			# if so, keep it!
			kept_signatures.append(most_enriched_signature)
			# next, exclude selected genes from further analysis: 1) adjust L 2) update set of excluded genes
			genes_used.update(most_enriched_signature.genes) # add selected genes to set of used genes
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
		print 'Filtering: kept %d / %d signatures.' %(len(kept_signatures),len(signatures)); sys.stdout.flush()

	return kept_signatures

def remove_redundant_signatures(new_signatures,previous_signatures,GO):
	if len(previous_signatures) == 0:
		return new_signatures
	kept_signatures = []
	previous_terms = set([enr.term[0] for enr in previous_signatures])
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

	mHG_X = go_mHG_X
	mHG_L = go_mHG_L
	if mHG_L == 0:
		mHG_L = len(genes)

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
	W = W[:,:test_pc] # truncate loading matrix
	final_signatures = []
	p = len(genes)
	res_var = None
	all_genes = set(genes)
	total_var = 0.0
	for pc in range(test_pc):

		print
		print '-'*70
		print "PC %d explains %.1f%% of the total variance." %(pc+1,100*frac[pc])
		total_var += frac[pc]
		print "The new cumulative fraction of total variance explained is %.1f%%." %(100*total_var)
		sys.stdout.flush()

		print "Testing for GO enrichment...", ; sys.stdout.flush()
		signatures_dsc = get_pc_signatures(M_enrich,W,pc+1,genes,mHG_X,mHG_L,go_pvalue_threshold,go_fold_enrichment_threshold)
		signatures_asc = get_pc_signatures(M_enrich,W,-pc-1,genes,mHG_X,mHG_L,go_pvalue_threshold,go_fold_enrichment_threshold)
		signatures = signatures_dsc + signatures_asc

		print "# signatures:",len(signatures); sys.stdout.flush()
		before = len(signatures)

		signatures = remove_redundant_signatures(signatures,final_signatures,GO)
		removed = before - len(signatures)
		print "Kept %d non-redundant signatures (%d removed)." %(len(signatures),removed)
	
		print_signatures(signatures,GO)
		final_signatures.extend(signatures)
		print "Total no. of signatures so far:", len(final_signatures); sys.stdout.flush()

		pc += 1

	print
	print '='*70
	print 'GO-PCA generated %d signatures:'
	print_signatures(final_signatures,GO)
	sys.stdout.flush()

	result = GOPCAResult(W=W,mHG_X=go_mHG_X,mHG_L=go_mHG_L,signatures=final_signatures)
	with open(args.output_file,'w') as ofh:
		pickle.dump(result,ofh,pickle.HIGHEST_PROTOCOL)

	print "Done!"
	return 0

if __name__ == '__main__':
	return_code = main(read_args_from_cmdline())
	sys.exit(return_code)
