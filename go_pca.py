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

import common

import sys
import argparse
import csv
import cPickle as pickle
import itertools as it

import numpy as np
from sklearn.decomposition import PCA

from tools import misc
import common
from go_enrichment import GOEnrichment

from goparser.parser import GOParser
from go_pca_objects import mHGTermResultWithPC,SignatureMatrix,NamedGeneSet

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-e','--expression-file',required=True)
	parser.add_argument('-g','--go-pickle-file',required=True)
	parser.add_argument('-x','--exclude-gene-file',default=None)

	parser.add_argument('-ag','--annotation-gene-file',required=True)
	parser.add_argument('-at','--annotation-term-file',required=True)
	parser.add_argument('-am','--annotation-matrix-file',required=True)

	parser.add_argument('-c','--max-components',type=int,default=0) # max. no. of PCA components to use
	parser.add_argument('-p','--p-value',type=float,default=0.01) # will be bonferroni-corrected
	parser.add_argument('-f','--fold-enrichment',type=float,default=1.5)
	parser.add_argument('-t','--top-genes',type=int,default=0) # by variance
	parser.add_argument('--min-pc-variance',type=float,default=1.0) # in %
	parser.add_argument('--pcs',type=int,default=0) # max. # of principal components to test

	parser.add_argument('--mHG-X',type=int,default=5)
	parser.add_argument('--mHG-L',type=int,default=1000)

	parser.add_argument('-o','--output-file',required=True)

	parser.add_argument('-j','--jobs',type=int,default=1)
	parser.add_argument('--verbose',action='store_true')

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

	enrichments_asc = go_enrichment_analysis(M,a,genes,X,L,p_value_threshold,fold_enrichment_threshold)
	enrichments_desc = go_enrichment_analysis(M,a[::-1],genes,X,L,p_value_threshold,fold_enrichment_threshold)

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

	min_pc_variance = args.min_pc_variance
	verbose = args.verbose

	# read expression data
	genes,samples,E = common.read_expression(args.expression_file)

	n,p = E.shape
	print "Expression matrix dimensions:",E.shape

	# remove excluded genes
	if args.exclude_gene_file is not None:
		excluded = np.zeros(n,dtype=np.bool_)
		excluded_genes = misc.read_single(args.exclude_gene_file)
		for g in excluded_genes:
			idx = misc.bisect_index(genes,g)
			excluded[idx] = True
		print 'Excluding %d genes!' %(np.sum(excluded)); sys.stdout.flush()
		sel = np.nonzero(np.invert(excluded))[0]
		genes = [genes[i] for i in sel]
		E = E[sel,:]
		n,p = E.shape

	# remove genes without gene expression
	sel = np.nonzero(np.amax(E,axis=1)>0)[0]
	genes = [genes[i] for i in sel]
	E = E[sel,:]
	
	# filter genes based on variance
	total_var = np.sum(np.var(E,axis=1,ddof=1))
	top_genes = args.top_genes
	n = E.shape[0]
	if top_genes > 0:
		full_var = total_var
		var = np.var(E,axis=1,ddof=1)
		a = np.argsort(var)[::-1]
		sel = np.zeros(n,dtype=np.bool_)
		sel[a[:top_genes]] = True
		sel = np.nonzero(sel)[0]
		E = E[sel,:]
		genes = [genes[i] for i in sel]
		total_var = np.sum(np.var(E,axis=1))
		print "Selected top %d genes based on variance (removing %.1f%% of overall variance)." \
				%(top_genes,100-100*(total_var/full_var)); sys.stdout.flush()

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
	annotation_files = [args.annotation_gene_file,args.annotation_term_file,args.annotation_matrix_file]
	M.read_annotations(*annotation_files)
	print "done!"; sys.stdout.flush()

	p_value_threshold = args.p_value / float(2*m)
	print "Bonferroni-corrected p-value threshold: %.1e" %(p_value_threshold); sys.stdout.flush()
	fold_enrichment_threshold = args.fold_enrichment
	final_enrichment = []
	N = len(genes)
	X = args.mHG_X
	L = args.mHG_L
	max_components = args.max_components
	res_var = None
	n,p = E.shape
	all_genes = set(genes)
	d = 0
	total_var = 0
	while True:
		if max_components > 0 and d >= max_components:
			break
		print "Performing PCA #%d..." %(d+1), ; sys.stdout.flush()
		P = PCA(n_components = 1)
		#P = PCA(n_components = 2)

		print "done!"; sys.stdout.flush()

		P.fit(E.T)
		C = P.components_
		frac = P.explained_variance_ratio_
		current_var = (1.0-total_var)*frac[0]
		if 100*current_var < min_pc_variance:
			print "Variance threshold crossed! Stopping..."; sys.stdout.flush()
			break

		total_var += current_var
		print 'Variance explained (total / current PC): %.1f%% / %.1f%%' %(100*total_var,100*current_var); sys.stdout.flush()

		enrichment = pc_go_enrichment(M,d+1,C[0,:],genes,X,L,p_value_threshold,fold_enrichment_threshold)

		print "# enriched terms:",len(enrichment); sys.stdout.flush()
		before = len(enrichment)

		enrichment = remove_redundant_terms(enrichment,final_enrichment,GO)
		removed = before - len(enrichment)
		print "Kept %d non-redundant enrichments (%d removed)." %(len(enrichment),removed)
	
		print_enrichment(enrichment,GO)
		final_enrichment.extend(enrichment)
		print "# total enriched terms:", len(final_enrichment); sys.stdout.flush()

		Y = P.transform(E.T)

		y = np.atleast_2d(Y[:,0]).T
		c = np.atleast_2d(C[0,:])

		E_remove = y.dot(c).T
		E = E - E_remove

		d += 1
		if args.pcs > 0 and d == args.pcs:
			print "Maximum number of PCs reached! Stopping..."; sys.stdout.flush()
			break

	print_enrichment(final_enrichment,GO)

	with open(args.output_file,'w') as ofh:
		pickle.dump(final_enrichment,ofh,pickle.HIGHEST_PROTOCOL)

	print "Done!"
	return 0

if __name__ == '__main__':
	return_code = main(read_args_from_cmdline())
	sys.exit(return_code)
