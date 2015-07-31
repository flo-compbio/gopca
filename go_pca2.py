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
	parser.add_argument('-f','--go-fold-enrichment-threshold',type=float,default=2.0)
	parser.add_argument('-Xf','--mHG-X-frac',type=float,default=0.25) # 0=off
	parser.add_argument('-Xm','--mHG-X-min',type=int,default=5) # 0=off
	parser.add_argument('-L','--mHG-L',type=int,default=1000) # 0=off
	parser.add_argument('--mfe-pvalue-threshold',type=float,default=1e-3)
	parser.add_argument('--ignore-mfe',action='store_true')

	# allow filtering to be disabled
	parser.add_argument('--disable-local-filter',action='store_true')
	parser.add_argument('--disable-global-filter',action='store_true')

	# variance filter
	parser.add_argument('-m','--most-variable-genes',type=int,default=0)
	
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

def print_signatures(signatures,GO,ignore_mfe=False):
	a = None
	maxlength = 40
	if not ignore_mfe:
		a = sorted(range(len(signatures)),key=lambda i: -signatures[i].mfe)
	else:
		a = sorted(range(len(signatures)),key=lambda i: -signatures[i].enrichment.fold_enrichment)

	for i in a:
		sig = signatures[i]
		#term = GO.terms[sig.term[0]]
		#goterm_genes = GO.get_goterm_genes(term.id)
		print sig.get_pretty_format(max_name_length=maxlength)

def get_pc_signatures(M,W,pc,genes,X_frac,X_min,L,pval_thresh,mfe_pval_thresh,filtering=True,fe_thresh=None,ignore_mfe=False,quiet=False):
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
	ranked_genes = [genes[i] for i in a]

	# find enriched GO terms using the XL-mHG test
	enriched_terms = M.test_enrichment(ranked_genes,pval_thresh,X_frac,X_min,L,quiet=quiet)
	if not enriched_terms:
		return []

	# calculate MFE for enriched GO terms
	q = len(enriched_terms)
	mfe = np.zeros(q,dtype=np.float64)
	k = np.zeros(q,dtype=np.int64)
	if not quiet:
		print 'Calculating maximum fold enrichment for each enriched term...', ; sys.stdout.flush()
	for i,enr in enumerate(enriched_terms):
		k[i], mfe[i] = enr.get_max_fold_enrichment(mfe_pval_thresh)
	if not quiet:
		print 'done!'; sys.stdout.flush()

	# filter enriched GO terms by strength of enrichment (if threshold is provided)
	if fe_thresh is not None:
		es = None
		if not ignore_mfe:
			es = mfe
		else:
			es = np.float64([enr.fold_enrichment for enr in enriched_terms])
		sel = np.nonzero(es >= fe_thresh)[0]
		enriched_terms = [enriched_terms[i] for i in sel]
		if not quiet:
			maximal = ' maximal'
			if ignore_mfe:
				maximal = ''
			print 'Enrichment filter: Kept %d / %d enriched terms with%s fold enrichment >= %.1fx.' %(sel.size,q,maximal,fe_thresh)

	# generate signatures
	signatures = []
	q = len(enriched_terms)
	for i,enr in enumerate(enriched_terms):
		sig_genes = set(enr.genes[:k[i]])
		#assert k[i] <= enr.K
		signatures.append(GOPCASignature(sig_genes,pc,mfe[i],enr))
	if not quiet:
		print 'Generated %d GO-PCA signatures based on the enriched GO terms.' %(q); sys.stdout.flush()

	# filter enriched GO terms (if filtering is enabled)
	if filtering:
		if not quiet:
			print 'Local filtering of signatures...', ; sys.stdout.flush()
		before = len(signatures)
		signatures = filter_signatures(signatures,M,ranked_genes,X_frac,X_min,L,\
				pval_thresh,mfe_pval_thresh,fe_thresh,ignore_mfe,quiet=quiet)
		if not quiet:
			print 'done!'
			print 'Local filter: kept %d / %d signatures.' %(len(signatures),before); sys.stdout.flush()

	return signatures

def filter_signatures(signatures,M,ranked_genes,X_frac,X_min,L,pval_thresh,mfe_pval_thresh,fe_thresh,ignore_mfe,quiet=False):

	if len(signatures) <= 1:
		return signatures

	# sort signatures by enrichment
	es = None
	q = len(signatures)
	if not ignore_mfe:
		es = np.float64([sig.mfe for sig in signatures])
	else:
		es = np.float64([sig.enrichment.fold_enrichment for sig in signatures])
	a = sorted(range(q), key=lambda i: -es[i])
	todo = [signatures[i] for i in a]
	print len(todo)

	# keep the most enriched signature
	most_enriched_sig = todo[0]
	kept_signatures = [most_enriched_sig]
	todo = todo[1:]

	# exclude genes in the most enriched signature
	ranked_genes = ranked_genes[:] # make a copy here!
	#genes_used = set(most_enriched_sig.genes)
	genes_used = set(most_enriched_sig.enrichment.genes)
	new_ranked_genes = []
	new_L = L
	for i,g in enumerate(ranked_genes):
		if g not in genes_used:
			new_ranked_genes.append(g)
		elif i < L: # gene was already used, adjust L if necessary
			new_L -= 1
	ranked_genes = new_ranked_genes
	L = new_L

	# start filtering
	K_max = max([sig.K for sig in todo])
	p = len(ranked_genes)
	mat = np.zeros((K_max+1,p+1),dtype=np.longdouble)
	while todo:
		most_enriched_sig = todo[0]
		term_id = most_enriched_sig.enrichment.term[0]

		# test if GO term is still enriched after removing all previously used genes
		enr = M.test_enrichment(ranked_genes,pval_thresh,X_frac,X_min,L,selected_term_ids=[term_id],mat=mat,quiet=True)
		assert len(enr) in [0,1]
		if enr: # enr will be an empty list if GO term does not meet the p-value threshold
			enr = enr[0]
			#print enr,'%d @ %d, s=%.1e' %(enr.mHG_k_n,enr.mHG_n,enr.mHG_s)

			# test fold enrichment threshold (if specified)
			still_enriched = False
			if fe_thresh is None:
				still_enriched = True
			elif not ignore_mfe:
				mfe = enr.get_max_fold_enrichment(mfe_pval_thresh)
				if mfe >= fe_thresh:
					still_enriched = True
			else:
				if enr.fold_enrichment >= fe_thresh:
					still_enriched = True

			if still_enriched:
				# if GO term is still considered enriched, keep it!
				kept_signatures.append(most_enriched_sig)
				# next, exclude selected genes from further analysis: 1) adjust L 2) update set of excluded genes
				before = len(genes_used)
				#genes_used.update(most_enriched_sig.genes) # add selected genes to set of used genes
				genes_used.update(most_enriched_sig.enrichment.genes) # add all genes of the GO term to set of used genes
				#print "test:",len(genes_used) - before
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

	return kept_signatures

def remove_redundant_signatures(new_signatures,previous_signatures,GO):
	if len(previous_signatures) == 0:
		return new_signatures
	kept_signatures = []
	previous_terms = set([sig.enrichment.term[0] for sig in previous_signatures])
	for sig in new_signatures:
		term_id = sig.enrichment.term[0]
		term = GO.terms[term_id]
		novel = True
		for t in set([term_id]) | term.ancestors | term.descendants:
			if t in previous_terms:
				novel = False
				break
		if novel:
			kept_signatures.append(sig)
	return kept_signatures

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	# input files
	expression_file = args.expression_file
	annotation_files = [args.annotation_gene_file,args.annotation_term_file,args.annotation_matrix_file]

	# parameters
	go_pvalue_threshold = args.go_pvalue_threshold
	go_fold_enrichment_threshold = args.go_fold_enrichment_threshold
	mHG_X_frac = args.mHG_X_frac
	mHG_X_min = args.mHG_X_min
	mHG_L = args.mHG_L

	mfe_pvalue_threshold = args.mfe_pvalue_threshold
	ignore_mfe = args.ignore_mfe

	most_variable_genes = args.most_variable_genes

	pc_num = args.pc_num
	pc_max = args.pc_max
	pc_cum_var = args.pc_cum_var/100.0

	disable_local_filter = args.disable_local_filter
	disable_global_filter = args.disable_global_filter

	# parameter checks?

	# misc options
	verbose = args.verbose

	# read expression data
	genes,samples,E = common.read_expression(expression_file)
	print "Expression matrix dimensions:", E.shape; sys.stdout.flush()

	# filter for most variable genes
	if most_variable_genes > 0:
		p = len(genes)
		sel = np.zeros(p,dtype=np.bool_)
		a = np.argsort(np.var(E,axis=1,ddof=1))
		a = a[::-1]
		sel[a[:most_variable_genes]] = True
		sel = np.nonzero(sel)[0]
		genes = [genes[i] for i in sel]
		E = E[sel,:]
		print 'Retained %d most variable genes.' %(most_variable_genes); sys.stdout.flush()
		print 'New expression matrix dimensions:', E.shape; sys.stdout.flush()

	if mHG_L == 0: # setting mHG_L to 0 will "turn off" the effect of the parameter (= set it to N)
		mHG_L = len(genes)

	# read GO data
	print "Loading GO term data...", ; sys.stdout.flush()
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

	# output cumulative fraction explained for each PC
	frac = M_pca.explained_variance_ratio_
	cum_frac = np.cumsum(frac)
	print 'Cumulative fraction of variance explained for first %d PCs:' %(compute_pc)
	print ', '.join(['%d: %.1f%%' %(pc+1,100*cum_frac[pc]) for pc in range(compute_pc)])
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

		#print "Testing for GO enrichment...", ; sys.stdout.flush()
		filtering = True
		if disable_local_filter:
			filtering = False
		signatures_dsc = get_pc_signatures(M_enrich,W,pc+1,genes,mHG_X_frac,mHG_X_min,mHG_L,go_pvalue_threshold,mfe_pvalue_threshold,\
				filtering,go_fold_enrichment_threshold,ignore_mfe=ignore_mfe)
		signatures_asc = get_pc_signatures(M_enrich,W,-pc-1,genes,mHG_X_frac,mHG_X_min,mHG_L,go_pvalue_threshold,mfe_pvalue_threshold,\
				filtering,go_fold_enrichment_threshold,ignore_mfe=ignore_mfe)
		signatures = signatures_dsc + signatures_asc

		print "# signatures:",len(signatures); sys.stdout.flush()
		before = len(signatures)

		if not disable_global_filter:
			signatures = remove_redundant_signatures(signatures,final_signatures,GO)
			print "Global filter: kept %d / %d signatures." %(len(signatures),before)
	
		print_signatures(signatures,GO,ignore_mfe)
		final_signatures.extend(signatures)
		print "Total no. of signatures so far:", len(final_signatures); sys.stdout.flush()

		pc += 1

	print
	print '='*70
	print 'GO-PCA generated %d signatures:' %(len(final_signatures))
	print_signatures(final_signatures,GO,ignore_mfe)
	sys.stdout.flush()

	result = GOPCAResult(genes=genes,W=W,mHG_X_frac=mHG_X_frac,mHG_X_min=mHG_X_min,mHG_L=mHG_L,signatures=final_signatures)
	with open(args.output_file,'w') as ofh:
		pickle.dump(result,ofh,pickle.HIGHEST_PROTOCOL)

	print "Done!"
	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
