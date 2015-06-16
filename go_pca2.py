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

	# additional input files
	parser.add_argument('--exclude-gene-file',default=None)

	# main parameters
	parser.add_argument('-t','--most-variable-genes',type=int,default=0) # filter for most variable genes (0=off)
	parser.add_argument('-p','--go-p-value-threshold',type=float,default=1e-6)
	parser.add_argument('-f','--go-fold-enrichment-threshold',type=float,default=1.5)

	parser.add_argument('-s','--seed',type=int,default=-1)

	# options for testing associations between genes and PCs
	parser.add_argument('-z','--gene-pc-zscore-file',default=None)
	parser.add_argument('--principal-components',type=int,default=0) # pre-determined (0=off)

	parser.add_argument('--pc-fdr',type=float,default=0.01) # only if z-score file is provided
	parser.add_argument('--pc-min-genes',type=int,default=100) # only if z-score file is provided (0=off)
	parser.add_argument('--pc-max',type=int,default=0) # hard limit (0=off)
	parser.add_argument('--pc-min-variance',type=float,default=0,help='in percent') # ignore PCs capturing less than this fraction of variance (0 = off)

	parser.add_argument('--mHG-X',type=int,default=5)
	parser.add_argument('--mHG-L',type=int,default=1000)

	#parser.add_argument('--disable-pc-test',action='store_true')
	#parser.add_argument('--pc-test-permutations',type=int,default=20)
	#parser.add_argument('--pc-test-permutation-percent',type=float,default=2.5)
	#parser.add_argument('--pc-test-jobs',type=int,default=1)

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

def pc_go_enrichment(M,pc,w,genes,X,L,p_value_threshold,fold_enrichment_threshold):
	a = np.argsort(w)

	L_asc = min(np.sum(w<0),L)
	print "Testing GO enrichment using ascending ordering of significant loadings (X=%d,L=%d)..." %(X,L_asc); sys.stdout.flush()
	enrichments_asc = go_enrichment_analysis(M,a,genes,X,L_asc,p_value_threshold,fold_enrichment_threshold)

	L_desc = min(np.sum(w>0),L)
	print "Testing GO enrichment using descending ordering of significant loadings (X=%d,L=%d)..." %(X,L_desc); sys.stdout.flush()
	enrichments_desc = go_enrichment_analysis(M,a[::-1],genes,X,L_desc,p_value_threshold,fold_enrichment_threshold)

	enrichments = []
	enrichments.extend(mHGTermResultWithPC.from_mHGTermResult(-pc,enr) for enr in enrichments_asc)
	enrichments.extend(mHGTermResultWithPC.from_mHGTermResult(pc,enr) for enr in enrichments_desc)

	return enrichments

def go_enrichment_analysis(M,ranking,genes,X,L,p_value_threshold,fold_enrichment_threshold):
	ranked_genes = [genes[i] for i in ranking]
	enrichment = M.test_enrichment(ranked_genes,X,L,quiet=True)

	enrichment = M.apply_thresholds(enrichment,p_value_threshold,fold_enrichment_threshold,quiet=False)
	enrichment = filter_enriched_go_terms(enrichment,M,ranking,genes,X,L,p_value_threshold,fold_enrichment_threshold,quiet=False)
	return enrichment

def filter_enriched_go_terms(enrichment,M,ranking,genes,X,L,p_value_threshold,fold_enrichment_threshold,quiet=True):

	todo = enrichment[:]
	genes_used = set()
	kept = []

	while todo:
		most_enriched = sorted(todo, key = lambda enr: -enr.fold_enrichment)[0]
		term_index = M.get_term_index(most_enriched.term[0])
		ranked_genes = [genes[i] for i in ranking if genes[i] not in genes_used]

		# test if still significant
		enr = M.test_enrichment(ranked_genes,X,L,selected_terms=[term_index],quiet=True)[0]
		if enr.p_value <= p_value_threshold and enr.fold_enrichment >= fold_enrichment_threshold:
			kept.append(most_enriched)
			genes_used.update(most_enriched.genes) # remove genes

		todo.remove(most_enriched)
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
	exclude_gene_file = args.exclude_gene_file
	zscore_file = args.gene_pc_zscore_file
	annotation_files = [args.annotation_gene_file,args.annotation_term_file,args.annotation_matrix_file]

	# parameters
	most_variable_genes = args.most_variable_genes
	go_p_value_threshold = args.go_p_value_threshold
	go_fold_enrichment_threshold = args.go_fold_enrichment_threshold
	pc_fdr = args.pc_fdr
	pc_min_genes = args.pc_min_genes
	pc_min_var = args.pc_min_variance/100.0
	mHG_X = args.mHG_X
	mHG_L = args.mHG_L

	principal_components = args.principal_components

	seed = args.seed
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
	modified = False

	if exclude_gene_file is not None:
		# remove excluded genes
		p = E.shape[0]
		sel = np.ones(p,dtype=np.bool_)
		exclude_genes = misc.read_single(exclude_gene_file)
		excluded = 0
		ignored = 0
		for g in exclude_genes:
			try:
				idx = misc.bisect_index(genes,g)
				sel[idx] = False
				excluded += 1
			except ValueError:
				ignored += 1
		if ignored > 0:
			print "Warning: Ignored %d unknown genes in excluded gene file." %(ignored); sys.stdout.flush()

		if excluded > 0:
			sel = np.nonzero(np.invert(excluded))[0]
			genes = [genes[i] for i in sel]
			E = E[sel,:]
			print 'Excluded %d genes!' %(excluded); sys.stdout.flush()
			modified = True

	# remove genes without gene expression
	sel = np.nonzero(np.amax(E,axis=1)>0)[0]
	not_expressed = len(genes) - sel.size
	if not_expressed > 0:
		genes = [genes[i] for i in sel]
		E = E[sel,:]
		print "Removed %d genes without expression." %(not_expressed); sys.stdout.flush()
		modified = True
	
	# filter genes based on variance
	total_var = np.sum(np.var(E,axis=1,ddof=1))
	p = E.shape[0]
	if most_variable_genes > 0:
		full_var = total_var
		var = np.var(E,axis=1,ddof=1)
		a = np.argsort(var)[::-1]
		sel = np.zeros(p,dtype=np.bool_)
		sel[a[:most_variable_genes]] = True
		sel = np.nonzero(sel)[0]
		E = E[sel,:]
		genes = [genes[i] for i in sel]
		total_var = np.sum(np.var(E,axis=1))
		print "Selected %d most variable genes (removing %.1f%% of the total variance)." \
				%(most_variable_genes,100-100*(total_var/full_var)); sys.stdout.flush()
		modified = True

	if modified:
		print "Final expression matrix has %d genes." %(len(genes))

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
	if zscore_file is not None:
		_,_,Z = common.read_gene_data(zscore_file)
		print E.shape[0],Z.shape[0]
		assert Z.shape[0] == E.shape[0] # make sure number of genes is correct

		# convert Z-scores to P-values (= one-sided tail)
		S = stats.norm.sf(Z,loc=0.0,scale=1.0)
		Z = None
		k = np.int64([fdr.fdr_bh(pval,pc_fdr)[0] for pval in S.T])
		print "TEST:",k

	# calculate principal components
	num_pcs = 50
	if principal_components > 0:
		num_pcs = principal_components
	elif S is not None:
		num_pcs = S.shape[1]
	P = RandomizedPCA(n_components = num_pcs)
	P.fit(E.T)
	frac = P.explained_variance_ratio_
	C = P.components_.T

	# run GO-PCA!
	final_enrichment = []
	p = len(genes)
	res_var = None
	all_genes = set(genes)
	total_var = 0.0
	d = C.shape[1]
	pc = 0
	while pc < d:

		#test_gene_pc_assoc(E,perc=)
		#break
		#current_var = (1.0-total_var)*frac[d]
		if frac[pc] < pc_min_var:
			print "Stopping! (PC does not explain enough % of total variance.)"; sys.stdout.flush()
			break

		c = C[:,pc]
		if S is not None:
			#pval = S[:,pc]
			k,crit_p = fdr.fdr_bh(S[:,pc],pc_fdr)
			if k < pc_min_genes: # test if sufficient genes are significantly associated with this PC
				print "Stopping! (PC has not enough genes significantly associated with it.)"; sys.stdout.flush()
				break

			# set non-significant loadings to zero
			not_significant = np.nonzero(S[:,pc] > crit_p)[0]
			#print "TEST:",crit_p,not_significant.size; sys.stdout.flush()
			c[not_significant] = 0.0

		s = ''
		if S is not None:
			s = 'has %d significant genes and ' %(k)
		print "PC %d %sexplains %.1f%% of the total variance." %(pc+1,s,100*frac[pc])
		total_var += frac[pc]
		print "The new cumulative fraction of total variance explained is %.1f%%." %(100*total_var)
		sys.stdout.flush()

		#pc += 1

		print "Testing for GO enrichment...", ; sys.stdout.flush()
		enrichment = pc_go_enrichment(M,pc+1,c,genes,mHG_X,mHG_L,go_p_value_threshold,go_fold_enrichment_threshold)

		print "# enriched terms:",len(enrichment); sys.stdout.flush()
		before = len(enrichment)

		enrichment = remove_redundant_terms(enrichment,final_enrichment,GO)
		removed = before - len(enrichment)
		print "Kept %d non-redundant enrichments (%d removed)." %(len(enrichment),removed)
	
		print_enrichment(enrichment,GO)
		final_enrichment.extend(enrichment)
		print "# total enriched terms:", len(final_enrichment); sys.stdout.flush()

		# manually remove PC
		#Y = P.transform(E.T)
		#y = np.atleast_2d(Y[:,0]).T
		#c = np.atleast_2d(C[0,:])
		#E_remove = y.dot(c).T
		#E = E - E_remove

		pc += 1

	print_enrichment(final_enrichment,GO)

	with open(args.output_file,'w') as ofh:
		pickle.dump(final_enrichment,ofh,pickle.HIGHEST_PROTOCOL)

	print "Done!"
	return 0

if __name__ == '__main__':
	return_code = main(read_args_from_cmdline())
	sys.exit(return_code)
