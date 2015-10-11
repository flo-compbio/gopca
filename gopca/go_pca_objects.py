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
import argparse
import re
import cPickle as pickle
import hashlib
import time
from copy import deepcopy
from collections import OrderedDict

import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import pearsonr

from gopca import common
from goparser import GOParser
from go_enrichment import GOEnrichment

class GOPCAArgumentParser(argparse.ArgumentParser):

	def __init__(self,*args,**kwargs):

		if 'description' not in kwargs or kwargs['description'] is None:
			kwargs['description'] = 'GO-PCA'
		argparse.ArgumentParser.__init__(self,*args,**kwargs)
		parser = self

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
		parser.add_argument('-D','--principal-components',type=int,default=0) # 0 = automatic

		###
		### Optional arguments
		###

		# log file
		parser.add_argument('-l','--log-file',default=None)

		# GO-PCA parameters
		parser.add_argument('-P','--pval-thresh',type=float,default=1e-6) # p-value threshold for GO enrichment
		parser.add_argument('-R','--sig-corr-thresh',type=float,default=0.5) # correlation threshold for signature genes
		parser.add_argument('-E','--escore-thresh',type=float,default=2.0) # XL-mHG enrichment score threshold

		parser.add_argument('-Xf','--mHG-X-frac',type=float,default=0.25) # 0=off
		parser.add_argument('-Xm','--mHG-X-min',type=int,default=5) # 0=off
		parser.add_argument('-L','--mHG-L',type=int,default=1000) # 0=off
		parser.add_argument('--escore-pval-thresh',type=float,default=1e-4) # p-value threshold for XL-mHG enrichment score calculation

		# variance filter
		parser.add_argument('-G','--select-variable-genes',type=int,default=0)

		# allow filtering to be disabled
		parser.add_argument('--disable-local-filter',action='store_true')
		parser.add_argument('--disable-global-filter',action='store_true')

		# output verbosity
		parser.add_argument('-v','--verbosity',type=int,default=3)

		# for automatically determining the number of PCs
		# (only used if -D is unset)
		parser.add_argument('-s','--seed',type=int,default=None)
		parser.add_argument('-pp','--pc-permutations',type=int,default=15)
		parser.add_argument('-pz','--pc-zscore-thresh',type=float,default=2.0)

		### legacy options
		parser.add_argument('--go-part-of-cc-only',action='store_true')


class GOPCAConfig(object):

	valid_params = set(['n_components','pval_thresh','sig_corr_thresh',\
			'mHG_X_frac','mHG_X_min','mHG_L',\
			'escore_pval_thresh','escore_thresh',\
			'disable_local_filter','disable_global_filter',\
			'seed','pc_permutations','pc_zscore_thresh',\
			'go_part_of_cc_only'])

	def __init__(self,logger,**kwargs):
		supplied_params = set(kwargs.keys())
		unknown_params = supplied_params - GOPCAConfig.valid_params
		for param in sorted(unknown_params):
			logger.warning('GO-PCA parameter "%s" is unknown and will be ignored.' %(param))

		# require that all parameters are provided
		for k in list(GOPCAConfig.valid_params):
			assert k in supplied_params

		# make sure the parameters are valid
		assert isinstance(kwargs['mHG_X_frac'],(int,float))
		assert 0.0 <= float(kwargs['mHG_X_frac']) <= 1.0
		assert isinstance(kwargs['mHG_X_min'],int)
		assert kwargs['mHG_X_min'] >= 0
		assert isinstance(kwargs['mHG_L'],int)
		assert kwargs['mHG_L'] >= 0
		assert isinstance(kwargs['pval_thresh'],(int,float))
		assert 0.0 < float(kwargs['pval_thresh']) <= 1.0
		assert isinstance(kwargs['escore_pval_thresh'],(int,float))
		assert 0.0 < float(kwargs['escore_pval_thresh']) <= 1.0
		assert isinstance(kwargs['escore_thresh'],(int,float))
		assert float(kwargs['escore_thresh']) >= 0.0
		assert isinstance(kwargs['seed'],int) and kwargs['seed'] >= 0
		assert isinstance(kwargs['pc_permutations'],int) and kwargs['pc_permutations'] > 0
		assert isinstance(kwargs['pc_zscore_thresh'],float)

		kwargs = dict([k,kwargs[k]] for k in list(GOPCAConfig.valid_params))
		self.__dict__.update(kwargs)

	@classmethod
	def read_config_file(cls):
		raise NotImplemented

	def write_config_file(self,output_file):
		raise NotImplemented

	def __repr__(self):
		return '<GOPCAConfig object (%s)>' %('; '.join(['%s=%s' %(k,getattr(self,str(k))) for k in sorted(self.valid_attrs)]))

	def __str__(self):
		return '<GOPCAConfig object with parameters: %s>' %(', '.join(['%s=%s' %(k,getattr(self,str(k))) for k in sorted(self.valid_attrs)]))

	def __hash__(self):
		return hash(repr(self))

	def __eq__(self):
		if type(self) is not type(other):
			return False
		if repr(self) == repr(other):
			return True
		else:
			return False

class GOPCA(object):

	def __init__(self,logger,config):

		assert isinstance(logger,common.Logger)
		assert isinstance(config,GOPCAConfig)
		self.logger = logger
		self.config = config

		#time_str = time.strftime('%Y-%m-%d_%H:%M:%S')
		#print 'Current time:',time_str

		self.genes = None
		self.samples = None
		self.E = None

		self.GO = None
		self.annotations = None

	@property
	def n_components(self):
		return self.config.n_components

	@property
	def D(self):
		return self.config.n_components

	@property
	def p(self):
		if self.E is None:
			return 0
		else:
			return E.shape[0]

	@property
	def n(self):
		if self.E is None:
			return 0
		else:
			return E.shape[1]

	@property
	def pval_thresh(self):
		return self.config.pval_thresh

	@property
	def sig_corr_thresh(self):
		return self.config.sig_corr_thresh

	@property
	def L(self):
		return self.config.mHG_L

	@property
	def X_frac(self):
		return self.config.mHG_X_frac

	@property
	def X_min(self):
		return self.config.mHG_X_min

	@property
	def escore_pval_thresh(self):
		return self.config.escore_pval_thresh

	@property
	def seed(self):
		return self.config.seed

	@property
	def pc_permutations(self):
		return self.config.pc_permutations

	@property
	def pc_zscore_thresh(self):
		return self.config.pc_zscore_thresh

	@property
	def escore_thresh(self):
		return self.config.escore_thresh

	def message(self,s,endline=True,flush=True):
		self.logger.message(s,endline,flush)

	def warning(self,s,endline=True,flush=True):
		self.logger.warning(s,endline,flush)

	def error(self,s,endline=True,flush=True):
		self.logger.error(s,endline,flush)

	def read_expression(self,expression_file):
		self.message('Reading expression...')

		# calculate hash
		hashval = hashlib.md5(open(expression_file,'rb').read()).hexdigest()
		self.message('MD5 hash: %s' %(hashval))

		genes,samples,E = common.read_expression(expression_file)
		p,n = E.shape
		self.message('Expression matrix size: p = %d genes x n = %d samples.' %(p,n))

		# adjust mHG L parameter, if necessary
		if self.config.mHG_L == 0 or self.config.mHG_L > p:
			self.config.mHG_L = p
			self.message('Set mHG_L to p (%d).'%(p))
		
		self.genes = genes
		self.samples = samples
		self.E = E

	def estimate_n_components(self,quiet=False):

		assert isinstance(self.E,np.ndarray) # make sure expression data has already been read

		# perform PCA
		if not quiet:
			self.message('Estimating the number of principal components (seed = %d)...' %(self.seed), endline=False)
		p,n = self.E.shape
		d_max = min(p,n-1)
		M_pca = PCA(n_components = d_max)
		M_pca.fit(self.E.T)

		d = M_pca.explained_variance_ratio_

		thresh = common.get_pc_explained_variance_threshold(self.E,self.pc_zscore_thresh,self.pc_permutations,self.seed)
		d_est = np.sum(d >= thresh)

		if not quiet:
			self.message('done!')
			self.message('The estimated number of PCs is %d.' %(d_est))
		self.config.n_components = d_est

	def print_signatures(self,signatures):
		a = None
		maxlength = 50
		a = sorted(range(len(signatures)),key=lambda i: -signatures[i].escore)

		for i in a:
			sig = signatures[i]
			self.message(sig.get_label(max_name_length=maxlength,include_pval=True))


	def filter_genes_by_variance(self,n_top):

		# find most variable genes
		p,n = self.E.shape
		sel = np.zeros(p,dtype=np.bool_)
		var = np.var(self.E,axis=1)
		a = np.argsort(var)
		a = a[::-1]
		sel[a[:n_top]] = True
		sel = np.nonzero(sel)[0]
		total_var = np.sum(var)

		# filtering
		self.genes = [self.genes[i] for i in sel]
		self.E = self.E[sel,:]

		# output some information
		lost_p = p - sel.size
		lost_var = total_var - np.sum(np.var(self.E,axis=1))
		self.message('Selected the %d most variable genes (excluded %.1f%% of genes, representing %.1f%% of total variance).' \
				%(n_top,100*(lost_p/float(p)),100*(lost_var/total_var)),flush=False)
		p,n = self.E.shape
		self.message('New expression matrix dimensions: %d genes x %d samples.' %(p,n))

		# adjust mHG L parameter, if necessary
		if self.config.mHG_L > sel.size:
			self.config.mHG_L = sel.size
			self.message('Adjusted "L" parameter to %d.' %(sel.size))

	def read_ontology(self,ontology_file):
		self.message('Reading ontology...')

		# calculate hash
		hashval = hashlib.md5(open(ontology_file,'rb').read()).hexdigest()
		self.message('(MD5 hash: %s)' %(hashval))

		self.GO = GOParser()
		self.GO.parse_ontology(ontology_file,part_of_cc_only=self.config.go_part_of_cc_only)

	def read_annotations(self,annotation_file):
		try:
			assert isinstance(self.E,np.ndarray)
		except AssertionError as e:
			self.error('Expression file not read yet!')
			raise e
		self.message('Reading annotations...',endline=False)

		# calculate hash
		hashval = hashlib.md5(open(annotation_file,'rb').read()).hexdigest()
		self.message('(MD5 hash: %s)' %(hashval),endline=False)

		self.annotations = common.read_annotations(annotation_file)

		n_assign = sum(len(v) for k,v in self.annotations.iteritems())
		self.message('(%d annotations) done!' %(n_assign))

	def run(self):
		assert isinstance(self.E,np.ndarray)
		assert isinstance(self.annotations,dict)
		if self.GO is not None:
			assert isinstance(self.GO,GOParser)

		t0 = time.time()

		# create GOEnrichment object
		self.message('Generating gene x GO term matrix...', endline=False)
		M_enrich = GOEnrichment(self.genes,self.annotations,self.logger)
		self.message('done!')

		# perform PCA
		self.message('Performing PCA...', endline=False)
		sys.stdout.flush()
		M_pca = PCA(n_components = self.n_components)
		Y = M_pca.fit_transform(self.E.T)
		self.message('done!')

		# output cumulative fraction explained for each PC
		frac = M_pca.explained_variance_ratio_
		cum_frac = np.cumsum(frac)
		self.message('Cumulative fraction of variance explained by the first %d PCs: %.1f%%' \
				%(self.n_components,100*cum_frac[-1]))

		# generate signatures
		W = M_pca.components_.T
		final_signatures = []
		p = len(self.genes)
		all_genes = set(self.genes)
		total_var = 0.0
		res_var = None
		for pc in range(self.n_components):

			total_var += frac[pc]
			self.message('',flush=False)
			self.message('-'*70,flush=False)
			self.message('PC %d explains %.1f%% of the variance.' %(pc+1,100*frac[pc]),flush=False)
			self.message('The new cumulative fraction of variance explained is %.1f%%.' %(100*total_var))

			#print "Testing for GO enrichment...", ; sys.stdout.flush()
			signatures_dsc = self.get_pc_signatures(M_enrich,W,pc+1)
			signatures_asc = self.get_pc_signatures(M_enrich,W,-pc-1)
			signatures = signatures_dsc + signatures_asc

			self.message("# signatures:",len(signatures))
			before = len(signatures)

			if not self.config.disable_global_filter:
				signatures = self.global_filter(signatures,final_signatures,self.GO)
				self.message("Global filter: kept %d / %d signatures." %(len(signatures),before))
		
			self.print_signatures(signatures)
			final_signatures.extend(signatures)
			self.message("Total no. of signatures so far: %d" %(len(final_signatures)))

			pc += 1

		self.message('')
		self.message('='*70)
		self.message('GO-PCA generated %d signatures:' %(len(final_signatures)))
		self.print_signatures(final_signatures)

		S = np.float64([common.get_signature_expression(self.genes,self.E,sig.genes) for sig in final_signatures])

		# it's important that we generate a copy of self.config for the result
		result = GOPCAResult(deepcopy(self.config),self.genes,self.samples,W,Y,final_signatures,S)

		t1 = time.time()
		self.message('GO-PCA runtime: %.2fs' %(t1-t0))

		return result

	def local_filter(self,M_enrich,enriched_terms,ranked_genes):
		# implements GO-PCA's "local" filter
		# returns enrichments that pass the filter

		if len(enriched_terms) <= 1:
			return enriched_terms

		# sort enriched terms by enrichment
		q = len(enriched_terms)
		a = sorted(range(q),key=lambda i:-enriched_terms[i].escore)
		todo = [enriched_terms[i] for i in a]

		# keep the most enriched term
		most_enriched = todo[0]
		kept_terms = [most_enriched]
		todo = todo[1:]

		# exclude genes annotated with the most enriched term
		genes_used = set(most_enriched.genes)
		new_ranked_genes = []
		L = self.L
		new_L = L
		for i,g in enumerate(ranked_genes):
			if g not in genes_used:
				new_ranked_genes.append(g)
			elif i < L: # gene was already used, adjust L if necessary
				new_L -= 1
		ranked_genes = new_ranked_genes
		L = new_L

		# start filtering
		K_max = max([enr.K for enr in todo])
		p = len(ranked_genes)
		mat = np.zeros((K_max+1,p+1),dtype=np.longdouble)
		while todo:
			most_enriched = todo[0]
			term_id = most_enriched.term[0]

			# test if GO term is still enriched after removing all previously used genes
			enr = M_enrich.get_enriched_terms(ranked_genes,self.pval_thresh,\
					self.X_frac,self.X_min,L,\
					escore_pval_thresh=self.escore_pval_thresh,selected_term_ids=[term_id],\
					mat=mat,verbosity=1)
			assert len(enr) in [0,1]
			if enr: # enr will be an empty list if GO term does not meet the p-value threshold
				enr = enr[0]
				#print enr,'%d @ %d, s=%.1e' %(enr.mHG_k_n,enr.mHG_n,enr.mHG_s)

				# test fold enrichment threshold (if specified)
				still_enriched = False
				if self.escore_thresh is None:
					still_enriched = True
				elif enr.escore >= self.escore_thresh:
					still_enriched = True

				if still_enriched:
					# if GO term is still considered enriched, keep it!
					kept_terms.append(most_enriched)
					# next, exclude selected genes from further analysis: 1) adjust L 2) update set of excluded genes
					genes_used.update(most_enriched.genes) # add all genes of the GO term to set of used genes
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

		self.message('Filtering: Kept %d / %d enriched terms.' \
				%(len(kept_terms),len(enriched_terms)),flush=True)
		return kept_terms

	def global_filter(self,new_signatures,previous_signatures,GO):
		if len(previous_signatures) == 0:
			return new_signatures
		kept_signatures = []
		previous_terms = set([sig.term[0] for sig in previous_signatures])
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

	def get_pc_signatures(self,M,W,pc):
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
		ranked_genes = [self.genes[i] for i in a]

		# find enriched GO terms using the XL-mHG test
		# get_enriched_terms also calculates the enrichment score, but does not use it for filtering
		enriched_terms = M.get_enriched_terms(ranked_genes,self.pval_thresh,self.X_frac,self.X_min,self.L,self.escore_pval_thresh)
		if not enriched_terms:
			return []

		# filter enriched GO terms by strength of enrichment (if threshold is provided)
		if self.escore_thresh is not None:
			q_before = len(enriched_terms)
			enriched_terms = [t for t in enriched_terms if t.escore >= self.escore_thresh]
			q = len(enriched_terms)
			self.message('Enrichment filter: Kept %d / %d enriched terms with enrichment score >= %.1fx.' \
					%(q,q_before,self.escore_thresh))

		# filter enriched GO terms (local filter)
		if not self.config.disable_local_filter:
			q_before = len(enriched_terms)
			enriched_terms = self.local_filter(M,enriched_terms,ranked_genes)
			q = len(enriched_terms)
			self.message('Local filter: Kept %d / %d enriched terms.' %(q,q_before))

		# generate signatures
		signatures = []
		q = len(enriched_terms)
		for j,enr in enumerate(enriched_terms):
			signatures.append(self.get_signature(pc,enr))
		self.message('Generated %d GO-PCA signatures based on the enriched GO terms.' %(q))

		return signatures

	def get_signature(self,pc,enr):
		""" Generate signature seed by selecting the X genes most correlated with average expression. """

		# select genes above cutoff giving rise to XL-mHG test statistic
		enr_genes = enr.genes[:enr.mHG_k_n]

		# calculate average expression
		indices = np.int64([self.genes.index(g) for g in enr_genes])
		E_enr = self.E[indices,:]
		mean = np.mean(common.get_standardized_matrix(E_enr),axis=0)

		# calculate seed based on the X genes most strongly correlated with the average
		corr = np.float64([pearsonr(mean,e)[0] for e in E_enr])
		a = np.argsort(corr)
		a = a[::-1]
		seed = np.mean(common.get_standardized_matrix(E_enr[a[:enr.X],:]),0)

		# select all other genes with correlation of at least sig_corr_thresh
		additional_indices = np.int64([i for i in a[enr.X:] if pearsonr(seed,E_enr[i,:])[0] > self.sig_corr_thresh])
		sel = np.r_[a[:enr.X],additional_indices]
		sig_genes = [enr_genes[i] for i in sel]
		sig_E = E_enr[sel,:]

		return GOPCASignature(sig_genes,sig_E,pc,enr)

class GOPCASignature(object):

	abbrev = [('positive ','pos. '),('negative ','neg. '),('interferon-','IFN-'),('proliferation','prolif.'),('signaling','signal.')]

	def __init__(self,genes,E,pc,xlmhg_result,label=None):
		self.genes = tuple(genes) # genes in the signature (NOT equal to self.enr.genes, which contains the gene names corresponding to all the 1's)
		self.E = E # expression of the genes in the signture, with ordering matching that of self.genes
		self.pc = pc # principal component (sign indicates whether ordering was ascending or descending)
		self.xlmhg_result = xlmhg_result # GO enrichment this signature is based on

	def __repr__(self):
		return '<GOPCASignature: label="%s", pc=%d, es=%.1f; %s>' \
				%(self.label,self.pc,self.escore,repr(self.enr))

	def __str__(self):
		return '<GO-PCA Signature "%s" (PC %d / E-score %.1fx / %s)>' \
				%(self.label,self.pc,self.escore, str(self.enr))

	def __hash__(self):
		return hash(repr(self))

	def __eq__(self,other):
		if type(self) is not type(other):
			return False
		elif repr(self) == repr(other):
			return True
		else:
			return False

	@property
	def enr(self):
		return self.xlmhg_result

	@property
	def term(self):
		""" The GO term that the signature is based on. """
		return self.enr.term

	@property
	def term_id(self):
		return self.enr.term[0]

	@property
	def term_name(self):
		return self.enr.term[3]

	@property
	def pval(self):
		""" The enrichment p-value of the GO term that the signature is based on. """
		return self.enr.pval

	@property
	def escore(self):
		return self.enr.escore
	
	@property
	def k(self):
		""" The number of genes in the signature. """
		return len(self.genes)

	@property
	def K(self):
		""" The number of genes annotated with the GO term whose enrichment led to the generation of the signature. """
		return self.enr.K

	@property
	def n(self):
		""" The threshold used to generate the signature. """
		return self.enr.ranks[self.k-1]

	@property
	def N(self):
		""" The total number of genes in the data. """
		return self.enr.N

	@property
	def label(self):
		return self.get_label(include_id=False)

	@property
	def median_correlation(self):
		C = np.corrcoef(self.E)
		ind = np.triu_indices(self.k,k=1)
		return np.median(C[ind])

	@property
	def escore_pval_thresh(self):
		return self.enr.escore_pval_thresh

	@property
	def gene_list(self):
		return ','.join(sorted(self.genes))

	def get_ordered_dict(self):
		elements = OrderedDict([
				['label',['Label',r'%s']],
				['pc',['PC',r'%d']],
				['term_id',['GO Term ID',r'%s']],
				['k',['k',r'%d']],['K',['K',r'%d']],
				['pval',['P-value',r'%.1e']],
				['escore',['E-score (psi=%.1e)' %(self.escore_pval_thresh),r'%.1f']],
				['median_correlation',['Median Correlation',r'%.2f']],
				['gene_list',['Genes',r'%s']]
		])
		od = OrderedDict([v[0],v[1] % (getattr(self,k))] for k,v in elements.iteritems())
		return od

	def get_label(self,max_name_length=0,include_stats=True,include_id=True,include_pval=False,include_collection=True):
		enr = self.enr

		term = enr.term
		term_name = term[3]
		for abb in self.abbrev:
			term_name = re.sub(abb[0],abb[1],term_name)
		if max_name_length > 0 and len(term_name) > max_name_length:
			term_name = term_name[:(max_name_length-3)] + '...'

		term_str = term_name
		if include_collection:
			term_str = '%s: %s' %(term[2],term_str)

		if include_id:
			term_str = term_str + ' (%s)' %(term[0])

		stats_str = ''
		if include_stats:
			if include_pval:
				stats_str = ' [%d:%d/%d, p=%.1e]' \
						%(self.pc,self.k,self.K,self.pval)
			else:
				stats_str = ' [%d:%d/%d]' \
						%(self.pc,self.k,self.K)

		return term_str + stats_str


class GOPCAResult(object):
	def __init__(self,config,genes,samples,W,Y,signatures,S):

		# W = PCA loading matrix
		# Y = PCA score matrix
		# S = GO-PCA signature matrix

		# checks
		assert isinstance(config,GOPCAConfig)
		assert isinstance(genes,list) or isinstance(genes,tuple)
		assert isinstance(samples,list) or isinstance(samples,tuple)
		assert isinstance(W,np.ndarray)
		assert isinstance(signatures,list) or isinstance(signatures,tuple)
		for s in signatures:
			assert isinstance(s,GOPCASignature)
		assert isinstance(S,np.ndarray)

		assert W.shape[0] == len(genes)
		assert Y.shape[0] == len(samples)
		assert W.shape[1] == Y.shape[1]
		assert S.shape[0] == len(signatures)
		assert S.shape[1] == len(samples)

		# initialization
		self.config = config
		self.genes = tuple(genes)
		self.samples = tuple(samples)
		self.W = W
		self.Y = Y
		self.signatures = tuple(signatures)
		self.S = S

	@property
	def p(self):
		return len(self.genes)

	@property
	def n(self):
		return len(self.samples)

	@property
	def D(self):
		return self.W.shape[1]

	@property
	def q(self):
		return len(self.signatures)

	@property
	def L(self):
		return self.config.mHG_L

	@property
	def X_frac(self):
		return self.config.mHG_X_frac

	@property
	def X_min(self):
		return self.config.mHG_X_min

	@property
	def escore_pval_thresh(self):
		return self.config.escore_pval_thresh

	@property
	def escore_thresh(self):
		return self.config.escore_thresh

	def save(self,output_file):
		#self.message('Saving GO-PCA result to file "%s"...' %(output_file),endline=False)
		with open(output_file,'wb') as ofh:
			pickle.dump(self,ofh,pickle.HIGHEST_PROTOCOL)
		#self.message("done!")

		#config = GOPCAConfig(**conf_dict)
		#result = GOPCAResult(config,genes,samples,W,final_signatures,S)
		#with open(output_file,'wb') as ofh:
		#	pickle.dump(result,ofh,pickle.HIGHEST_PROTOCOL)

	@staticmethod
	def load(file_name):
		result = None
		with open(file_name,'rb') as fh:
			result = pickle.load(fh)
		return result

	def __repr__(self):
		conf_hash = hash(self.config)
		gene_hash = hash(self.genes)
		sample_hash = hash(self.samples)
		sig_hash = hash((hash(sig) for sig in self.signatures))
		return '<GOPCAResult object (config hash: %d; gene hash: %d; sample hash: %d; # PCs: %d; signatures: %d; signature hash: %d)>' \
				%(conf_hash,gene_hash,sample_hash,self.D,self.q,sig_hash)

	def __str__(self):
		conf = self.config
		return '<GOPCAResult object (%d signatures); mHG parameters: X_frac=%.2f, X_min=%d, L=%d; \
				# genes (p) = %d, # principal components (D) = %d>' \
				%(self.q,conf.mHG_X_frac,conf.mHG_X_min,cof.mHG_L,self.p,self.D)

	def __eq__(self,other):
		if type(self) is not type(other):
			return False
		elif repr(self) == repr(other):
			return True
		else:
			return False


