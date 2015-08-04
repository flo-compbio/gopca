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

#from go_enrichment import mHGTermResult

class GOPCAResult(object):
	def __init__(self,genes,W,mHG_X_frac,mHG_X_min,mHG_L,signatures):
		assert len(genes) == W.shape[0]
		self.genes = genes
		self.W = W
		self.mHG_X_frac = mHG_X_frac
		self.mHG_X_min = mHG_X_min
		self.mHG_L = mHG_L
		self.signatures = signatures

	def __repr__(self):
		sig_hash = hash((hash(sig) for sig in self.signatures))
		p,d = self.W.shape
		return "<GO-PCA result (mHG_X_frac=%.2f, mHG_X_min=%d, mHG_L=%d); loading matrix dimensions: (%d,%d); hash of signature list: %d>" \
				%(self.mHG_X_frac,self.mHG_X_min,self.mHG_L,p,d,sig_hash)

	def __str__(self):
		q = len(self.signatures)
		p,d = self.W.shape
		return "<GO-PCA result (%d signatures); mHG parameters: X_frac=%.2f, X_min=%d, L=%d; # genes (p) = %d, # principal components (d) = %d>" \
				%(q,self.mHG_X_frac,self.mHG_X_min,self.mHG_L,p,d)

	def __eq__(self,other):
		if type(self) != type(other):
			return False
		elif repr(self) == repr(other):
			return True
		else:
			return False

class GOPCASignature(object):

	abbrev = [('positive ','pos. '),('negative ','neg. '),('interferon-','IFN-'),('proliferation','prolif.'),('signaling','signal.')]

	def __init__(self,genes,pc,mfe,enrichment,label=None):
		self.genes = set(genes) # genes in the signature
		self.pc = pc # principal component (sign indicates whether ordering was ascending or descending)
		self.mfe = mfe # maximum fold enrichment
		self.enrichment = enrichment # GO enrichment this signature is based on
		if label is None:
			enr = enrichment
			label = '%s: %s (%d:%d/%d)' %(enr.term[2],enr.term[3],pc,enr.k,enr.K)
		self.label = label # signature label

	def __repr__(self):
		return '<GOPCASignature: label="%s", pc=%d, mfe=%.1f; %s>' \
				%(self.label,self.pc,self.mfe,repr(self.enrichment))

	def __str__(self):
		return '<GO-PCA Signature "%s" (PC %d / MFE %.1fx / %s)>' \
				%(self.label,self.pc,self.mfe, str(self.enrichment))

	def __hash__(self):
		return hash(repr(self))

	def __eq__(self,other):
		if type(self) != type(other):
			return False
		elif repr(self) == repr(other):
			return True
		else:
			return False

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

	def get_pretty_format(self,omit_acc=False,nitty_gritty=True,max_name_length=0):
		enr = self.enrichment

		term = enr.term
		term_name = term[3]
		for abb in self.abbrev:
			term_name = re.sub(abb[0],abb[1],term_name)
		if len(term_name) > max_name_length:
			term_name = term_name[:(max_name_length-3)] + '...'

		term_str = '%s: %s' %(term[2],term_name)
		if not omit_acc:
			term_str = term_str + ' (%s)' %(term[0])

		details = ''
		if nitty_gritty:
			details = ' [%d/%d genes,pc=%d,pval=%.1e,mfe=%.1fx]' \
					%(self.k,self.K,self.pc,self.pval,self.mfe)

		return '%s%s' %(term_str,details)

	def get_pretty_format_GO(self,GO,omit_acc=False,nitty_gritty=True,max_name_length=0):
		enr = self.enrichment
		term = GO.terms[enr.term[0]]
		goterm_genes = GO.get_goterm_genes(term.id)
		details = ''
		if nitty_gritty:
			details = ' [%d/%d genes,n=%d,pc=%d,mfe=%.1fx,pval=%.1e]' \
					%(self.k,self.K,self.n,self.pc,self.mfe,self.pval)
		return '%s%s' %(term.get_pretty_format(omit_acc=omit_acc,max_name_length=max_name_length),details)

	#@staticmethod
	#def from_mHGTermResult(result,pc,es=None):
	#	return GOPCASignature(result.term,result.p_value,result.N,result.n,result.K,result.genes,pc,es)

c="""
class GeneSet(object):
	def __init__(self,genes):
		self.genes = frozenset(genes)

	def __repr__(self):
		return '<GeneSet with %d genes (hash:%d)' %(len(self.genes),hash(self.genes))

	def __str__(self):
		return '<GeneSet with %d genes: %s>' %(len(self.genes),', '.join(sorted(self.genes)))

	def __eq__(self,other):
		if type(self) != type(other):
			return False
		elif repr(self) == repr(other):
			return True
		else:
			return False


class NamedGeneSet(GeneSet):
	def __init__(self,name,genes):
		GeneSet.__init__(self,genes)
		self.name = name

	def __repr__(self):
		return '<NamedGeneSet "%s" with %d genes (hash:%d)' %(self.name,len(self.genes),hash(self.genes))

	def __str__(self):
		return '<NamedGeneSet "%s" with %d genes: %s>' %(self.name,len(self.genes),', '.join(sorted(self.genes)))

	def __eq__(self,other):
		if type(self) != type(other):
			return False
		elif repr(self) == repr(other):
			return True
		else:
			return False

	@staticmethod
	def from_GeneSet(name,geneset):
		return NamedGeneSet(name,geneset.genes)

class SignatureMatrix(object):
	def __init__(self,genesets,S):
		assert len(genesets) == S.shape[0]
		self.S = S
		self.genesets = genesets

	def save(self,fn):
		with open(fn,'w') as ofh:
			pickle.dump(self,ofh,pickle.HIGHEST_PROTOCOL)
"""

