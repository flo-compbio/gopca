#!/usr/bin/env python2.7

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

import argparse
import csv
import cPickle as pickle
import itertools as it
import hashlib
import time

import numpy as np
from scipy import stats
from sklearn.decomposition import PCA

from genometools import misc
from goparser import GOParser
from gopca.common import Logger,print_signatures
from gopca.go_enrichment import GOEnrichment
from gopca.go_pca_objects import GOPCAConfig,GOPCASignature,GOPCAResult,GOPCA
from gopca.printf import printf

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='GO-PCA')

	###
	### Required arguments
	###

	# input files
	parser.add_argument('-e','--expression-file',required=True)
	parser.add_argument('-a','--annotation-file',required=True)
	parser.add_argument('-t','--ontology-file',default=None)

	# output file
	parser.add_argument('-o','--output-file',required=True)

	# number of principal components to test
	parser.add_argument('-D','--principal-components',type=int,required=True)

	###
	### Optional arguments
	###

	# log file
	parser.add_argument('-l','--log-file',default=None)

	# GO-PCA parameters
	parser.add_argument('-p','--pval-thresh',type=float,default=1e-6) # p-value threshold for GO enrichment
	parser.add_argument('-f','--msfe-thresh',type=float,default=2.0) # MSFE (fold enrichment) threshold
	parser.add_argument('-r','--sig-corr-threshold',type=float,default=0.5) # correlation threshold for signature genes

	parser.add_argument('-Xf','--mHG-X-frac',type=float,default=0.25) # 0=off
	parser.add_argument('-Xm','--mHG-X-min',type=int,default=5) # 0=off
	parser.add_argument('-L','--mHG-L',type=int,default=1000) # 0=off
	parser.add_argument('--msfe-pval-thresh',type=float,default=1e-4) # p-value threshold for MSFE calculation

	# variance filter
	parser.add_argument('-n','--select-variable-genes',type=int,default=0)

	# allow filtering to be disabled
	parser.add_argument('--disable-local-filter',action='store_true')
	parser.add_argument('--disable-global-filter',action='store_true')

	### legacy options
	parser.add_argument('--go-part-of-cc-only',action='store_true')

	# output verbosity
	#parser.add_argument('-v','--verbose',action='store_true')
	parser.add_argument('-v','--verbose',action='store_false')

	return parser.parse_args()

class GOPCA(object):

	def __init__(logger,config,ontology_file,log_file,verbosity=3):

		self.config = config
		self.logger = common.Logger(verbosity,log_file = log_file)

		#time_str = time.strftime('%Y-%m-%d_%H:%M:%S')
		#print 'Current time:',time_str


		self.genes = None
		self.samples = None
		self.E = None

		self.GO = None
		self.annotations = None

		self.result = None

	def message(self,s,endline=True,flush=True):
		self.logger.message(s,endline,flush)

	def warning(self,s,endline=True,flush=True):
		self.logger.warning(s,endline,flush)

	def error(self,s,endline=True,flush=True):
		self.logger.error(s,endline,flush)

	def read_expression(self,expression_file):
		self.message('Reading expression...',endline=False)

		expression_hash = hashlib.md5(open(expression_file,'rb').read()).hexdigest()
		print 'Expression file MD5 hash: %s' %(expression_hash)

		genes,samples,E = common.read_expression(expression_file)
		self.genes = genes
		self.samples = samples
		self.E = E

		self.message('done!')

	def filter_genes_by_variance(self,n_top):
		genes = self.genes
		E = self.E

		# find most variable genes
		p,n = E.shape
		sel = np.zeros(p,dtype=np.bool_)
		var = np.var(E,axis=1)
		a = np.argsort(var)
		a = a[::-1]
		sel[a[:n_top]] = True
		sel = np.nonzero(sel)[0]
		total_var = np.sum(var)

		# filtering
		self.genes = [genes[i] for i in sel]
		self.E = E[sel,:]

		# output some information
		lost_p = p - sel.size
		lost_var = total_var - np.sum(np.var(self.E,axis=1))
		self.message('Retained the %d most variable genes (excluded %.1f%% of genes, representing %.1f%% of total variance).' \
				%(n_top,100*(lost_p/float(p)),100*(lost_var/total_var)),flush=False)
		self.message('New expression matrix dimensions:i ' + str(self.E.shape))

	def read_ontology(self,gene_file,ontology_file):
		self.message('Reading ontology...',endline=False)

		# hash

		self.GO = GOParser()
		self.GO.parse_ontology(ontology_file,part_of_cc_only=self.go_part_of_cc_only)

		self.message('done!')

	def read_annotations(self,annotation_file):
		self.message('Reading annotations...',endline=False)

		# hash

		self.annotations = common.read_annotations(annotation_file)

		self.message('done!')

def get_pc_signatures(M,W,pc,E,genes,X_frac,X_min,L,pval_thresh,msfe_pval_thresh=None,filtering=True,msfe_thresh=None,quiet=False):
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

	# calculate MSFE for enriched GO terms
	q = len(enriched_terms)
	msfe = np.zeros(q,dtype=np.float64)
	if not quiet:
		print 'Calculating maximum fold enrichment for each enriched term...', ; sys.stdout.flush()
	for i,enr in enumerate(enriched_terms):
		_, msfe[i] = enr.get_max_sig_fold_enrichment(msfe_pval_thresh)
	if not quiet:
		print 'done!'; sys.stdout.flush()

	# filter enriched GO terms by strength of enrichment (if threshold is provided)
	if msfe_thresh is not None:
		sel = np.nonzero(msfe >= msfe_thresh)[0]
		enriched_terms = [enriched_terms[i] for i in sel]
		if not quiet:
			print 'Enrichment filter: Kept %d / %d enriched terms with maximal fold enrichment >= %.1fx.' %(sel.size,q,msfe_thresh)

	# filter enriched GO terms
	q_before = len(enriched_terms)
	enriched_terms = M.filter_enriched_terms(enriched_terms,ranked_genes,pval_thresh,X_frac,X_min,L,msfe_pval_thresh,msfe_thresh)
	q = len(enriched_terms)
	print 'Local filter: Kept %d / %d enriched terms.' %(q,q_before); sys.stdout.flush()

	signatures = []
	q = len(enriched_terms)
	for j,enr in enumerate(enriched_terms):
		sig_genes = enr.genes[:enr.mHG_k_n] # important!
		indices = np.zeros(len(sig_genes),dtype=np.int64)
		for i,g in enumerate(sig_genes):
			indices[i] = genes.index(g)
		sig_E = E[indices,:]
		signatures.append(GOPCASignature(sig_genes,sig_E,pc,msfe[j],enr))
	if not quiet:
		print 'Generated %d GO-PCA signatures based on the enriched GO terms.' %(q); sys.stdout.flush()

	# filter enriched GO terms (if filtering is enabled)
	if filtering:
		if not quiet:
			print 'Local filtering of signatures...', ; sys.stdout.flush()
		before = len(signatures)
		signatures = filter_signatures(signatures,M,ranked_genes,X_frac,X_min,L,\
				pval_thresh,msfe_pval_thresh,msfe_thresh,quiet=quiet)
		if not quiet:
			print 'done!'
			print 'Local filter: kept %d / %d signatures.' %(len(signatures),before); sys.stdout.flush()

	return signatures

def filter_signatures(signatures,M,ranked_genes,X_frac,X_min,L,pval_thresh,msfe_pval_thresh,msfe_thresh,quiet=False):

	if len(signatures) <= 1:
		return signatures

	# sort signatures by enrichment
	es = None
	q = len(signatures)
	es = np.float64([sig.msfe for sig in signatures])
	a = sorted(range(q), key=lambda i: -es[i])
	todo = [signatures[i] for i in a]
	#print len(todo)

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
			if msfe_thresh is None:
				still_enriched = True
			else:
				msfe = enr.get_max_sig_fold_enrichment(msfe_pval_thresh)
				if msfe >= msfe_thresh:
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
	ontology_file = args.ontology_file
	annotation_file = args.annotation_file

	# output file
	output_file = args.output_file

	# log file
	log_file = args.log_file

	# GO-PCA parameters
	D = args.principal_components
	pval_thresh = args.pval_thresh
	msfe_thresh = args.msfe_thresh
	mHG_X_frac = args.mHG_X_frac
	mHG_X_min = args.mHG_X_min
	mHG_L = args.mHG_L
	msfe_pval_thresh = args.msfe_pval_thresh

	sel_var_genes = args.select_variable_genes

	disable_local_filter = args.disable_local_filter
	disable_global_filter = args.disable_global_filter

	part_of_cc_only = args.go_part_of_cc_only
	verbose = args.verbose

	# make sure input files exist
	assert os.path.isfile(expression_file)
	assert os.path.isfile(annotation_file)
	if ontology_file is not None:
		assert os.path.isfile(ontology_file)

	# disable global filter if no ontology is provided
	if ontology_file is None:
		print 'Warning: Disabling global filter, since not ontology file was provided.'; sys.stdout.flush()
		disable_global_filter = True


	# generate GOPCAConfig object
	conf_params = ['D','mHG_X_frac','mHG_X_min','mHG_L',\
			'pval_thresh','msfe_pval_thresh','msfe_thresh',\
			'disable_local_filter','disable_global_filter'\
			'go_part_of_cc_only']
	conf_dict = dict([[k,locals()[k]] for k in conf_params])
	C = GOPCAConfig(conf_dict)
	M = GOPCA(config=C,log_file=log_file)

	# read expression data
	M.read_expression(expression_file)

	# read ontology
	if ontology_file is not None:
		M.read_ontology(ontology_file)

	# read annotations
	M.read_annotations(annotation_file)

	# read expression data
	print 'Reading expression data...', ; sys.stdout.flush()
	genes,samples,E = common.read_expression(expression_file)
	print 'done!'; sys.stdout.flush()
	print "Expression matrix dimensions:", E.shape; sys.stdout.flush()

	print 'Reading annotations...', ; sys.stdout.flush()
	annotation_hash = hashlib.md5(open(annotation_file,'rb').read()).hexdigest()
	annotations = common.read_annotations(annotation_file)
	n_assign = sum(len(v) for k,v in annotations.iteritems())
	print 'done!'; sys.stdout.flush()
	print 'Annotation file MD5 hash: %s' %(annotation_hash)
	print 'Read %d annotations for %d GO terms.' %(n_assign,len(annotations))

	# parse ontology
	if ontology_file is not None:
		print 'Reading ontology...', ; sys.stdout.flush()
		ontology_hash = hashlib.md5(open(ontology_file,'rb').read()).hexdigest()
		GO = GOParser()
		GO.parse_ontology(ontology_file,part_of_cc_only=part_of_cc_only,quiet=True)
		print 'done!'; sys.stdout.flush()
		print 'Ontology file MD5 hash: %s' %(ontology_hash)
		print 'Read ontology with %d terms.' %(len(GO.terms))

	# filter for most variable genes
	if var_filter_genes > 0:
	if mHG_L == 0: # setting mHG_L to 0 will "turn off" the effect of the parameter (= set it to N)
		mHG_L = len(genes)

	# sort genes alphabetically
	#gene_order = np.int64(misc.argsort(genes))
	#genes = [genes[i] for i in gene_order]
	#E = E[gene_order,:]

	# create GOEnrichment object
	print 'Generating gene x GO term matrix...', ; sys.stdout.flush()
	M_enrich = GOEnrichment(genes,annotations)
	print 'done!'; sys.stdout.flush()

	# perform PCA
	print 'Performing PCA...', ; sys.stdout.flush()
	sys.stdout.flush()
	M_pca = PCA(n_components = pc_num)
	M_pca.fit(E.T)
	print 'done!'

	# output cumulative fraction explained for each PC
	frac = M_pca.explained_variance_ratio_
	cum_frac = np.cumsum(frac)
	print 'Cumulative fraction of variance explained by the first %d PCs: %.1f%%' %(pc_num,100*cum_frac[-1])
	#print ', '.join(['%d: %.1f%%' %(pc+1,100*cum_frac[pc]) for pc in range(compute_pc)])
	#sys.stdout.flush()

	# run GO-PCA!
	W = M_pca.components_.T
	#W = W[:,:pc_num] # truncate loading matrix
	final_signatures = []
	p = len(genes)
	res_var = None
	all_genes = set(genes)
	total_var = 0.0
	filtering = not disable_local_filter
	for pc in range(pc_num):

		print
		print '-'*70
		print "PC %d explains %.1f%% of the variance." %(pc+1,100*frac[pc])
		total_var += frac[pc]
		print "The new cumulative fraction of variance explained is %.1f%%." %(100*total_var)
		sys.stdout.flush()

		#print "Testing for GO enrichment...", ; sys.stdout.flush()
		signatures_dsc = get_pc_signatures(M_enrich,W,pc+1,E,genes,mHG_X_frac,mHG_X_min,mHG_L,pval_thresh,msfe_pval_thresh,\
				filtering,msfe_thresh)
		signatures_asc = get_pc_signatures(M_enrich,W,-pc-1,E,genes,mHG_X_frac,mHG_X_min,mHG_L,pval_thresh,msfe_pval_thresh,\
				filtering,msfe_thresh)
		signatures = signatures_dsc + signatures_asc

		print "# signatures:",len(signatures); sys.stdout.flush()
		before = len(signatures)

		if not disable_global_filter:
			signatures = remove_redundant_signatures(signatures,final_signatures,GO)
			print "Global filter: kept %d / %d signatures." %(len(signatures),before)
	
		print_signatures(signatures,GO)
		final_signatures.extend(signatures)
		print "Total no. of signatures so far:", len(final_signatures); sys.stdout.flush()

		pc += 1

	print
	print '='*70
	print 'GO-PCA generated %d signatures:' %(len(final_signatures))
	print_signatures(final_signatures,GO)
	sys.stdout.flush()

	print 'Writing results to file "%s"...' %(output_file), ; sys.stdout.flush()
	config = GOPCAConfig(**conf_dict)
	S = np.float64([common.get_signature_expression(genes,E,sig.genes) for sig in final_signatures])
	result = GOPCAResult(config,genes,samples,W,final_signatures,S)
	with open(output_file,'wb') as ofh:
		pickle.dump(result,ofh,pickle.HIGHEST_PROTOCOL)
	print "done!"; sys.stdout.flush()

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
