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

import re
import cPickle as pickle

import numpy as np

from gopca.common import Logger

		
class GOPCAConfig(object):

	valid_params = set(['D',\
			'mHG_X_frac','mHG_X_min','mHG_L',\
			'pval_thresh','msfe_pval_thresh','msfe_thresh',\
			'disable_local_filter','disable_global_filter',\
			'go_part_of_cc_only'])

	def __init__(self,**kwargs):
		supplied_params = set(kwargs.keys())
		unknown_params = supplied_params - self.valid_params
		for param in sorted(unknown_params):
			print 'Warning: GO-PCA parameter "%s" is unknown (will be ignored).' %(k)

		# require that all parameters are specified
		for k in list(self.valid_params):
			assert k in supplied_params

		# to-do: make sure parmaeters parameter are valid
		kwargs = dict([k,kwargs[k]] for k in list(self.valid_params))
		self.__dict__.update(kwargs)


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


class GOPCASignature(object):

	abbrev = [('positive ','pos. '),('negative ','neg. '),('interferon-','IFN-'),('proliferation','prolif.'),('signaling','signal.')]

	def __init__(self,genes,E,pc,msfe,enrichment,label=None):
		self.genes = tuple(genes) # genes in the signature (NOT equal to self.enrichment.genes, which contains the gene names corresponding to all the 1's)
		self.E = E # expression of the genes in the signture, with ordering matching that of self.genes
		self.pc = pc # principal component (sign indicates whether ordering was ascending or descending)
		self.msfe = msfe # maximum significant fold enrichment
		self.enrichment = enrichment # GO enrichment this signature is based on
		if label is None:
			enr = enrichment
			label = '%s: %s (%d:%d/%d)' %(enr.term[2],enr.term[3],pc,enr.k,enr.K)
		self.label = label # signature label

	def __repr__(self):
		return '<GOPCASignature: label="%s", pc=%d, msfe=%.1f; %s>' \
				%(self.label,self.pc,self.msfe,repr(self.enrichment))

	def __str__(self):
		return '<GO-PCA Signature "%s" (PC %d / MSFE %.1fx / %s)>' \
				%(self.label,self.pc,self.msfe, str(self.enrichment))

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
	def term(self):
		return self.enrichment.term

	@property
	def pval(self):
		""" The enrichment p-value of the GO term that the signature is based on. """
		return self.enrichment.pval
	
	@property
	def k(self):
		""" The number of genes in the signature. """
		return len(self.genes)

	@property
	def K(self):
		""" The number of genes annotated with the GO term whose enrichment led to the generation of the signature. """
		return self.enrichment.K

	@property
	def n(self):
		""" The threshold used to generate the signature. """
		return self.enrichment.ranks[self.k-1]

	@property
	def N(self):
		""" The total number of genes in the data. """
		return self.enrichment.N

	def get_label(self,max_name_length=0,include_stats=True,include_id=True,include_pval=False,include_collection=True):
		enr = self.enrichment

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
	def __init__(self,config,genes,samples,W,signatures,S):

		# W = loading matrix
		# S = signature matrix

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
		assert S.shape[0] == len(signatures)
		assert S.shape[1] == len(samples)

		# initialization
		self.gopca_input = gopca_input
		self.genes = tuple(genes)
		self.samples = tuple(samples)
		self.W = W
		self.signatures = tuple(signatures)
		self.S = S

	@property
	def p(self):
		return len(self.genes)

	@property
	def n(self):
		return len(self.samples)

	@property
	def d(self):
		return self.W.shape[1]

	@property
	def q(self):
		return len(self.signatures)
		
	def __repr__(self):
		conf_hash = hash(self.config)
		gene_hash = hash(self.genes)
		sample_hash = hash(self.samples)
		sig_hash = hash((hash(sig) for sig in self.signatures))
		return '<GOPCAResult object (config hash: %d; gene hash: %d; sample hash: %d; # PCs: %d; signatures: %d; signature hash: %d)>' \
				%(conf_hash,gene_hash,sample_hash,self.d,self.q,sig_hash)

	def __str__(self):
		conf = self.config
		return '<GOPCAResult object (%d signatures); mHG parameters: X_frac=%.2f, X_min=%d, L=%d; \
				# genes (p) = %d, # principal components (d) = %d>' \
				%(self.q,conf.mHG_X_frac,conf.mHG_X_min,cof.mHG_L,self.p,self.d)

	def __eq__(self,other):
		if type(self) is not type(other):
			return False
		elif repr(self) == repr(other):
			return True
		else:
			return False
