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

import cPickle as pickle

from go_enrichment import mHGTermResult

class GOPCAResult(object):
	def __init__(self,genes,W,mHG_X,mHG_L,signatures):
		self.genes = genes
		self.W = W
		self.mHG_X = mHG_X
		self.mHG_L = mHG_L
		self.signatures = signatures

	def __repr__(self):
		sig_hash = hash((hash(sig) for sig in self.signatures))
		p,d = self.W.shape
		return "<GO-PCA result (mHG_X=%d, mHG_L=%d); loading matrix dimensions: (%d,%d); hash of signature list: %d>" %(mHG_X,mHG_L,p,d,sig_hash)

	def __str__(self):
		q = len(self.signatures)
		p,d = self.W.shape
		return "<GO-PCA result; mHG parameters: X=%d, L=%d; # genes (p) = %d, # principal components (d) = %d; %d signatures>" \
				%(mHG_X,mHG_L,p,d,q)

	def __eq__(self,other):
		if type(self) != type(other):
			return False
		elif repr(self) == repr(other):
			return True
		else:
			return False

class GOPCASignature(mHGTermResult):
	"""
	Stores GO-PCA signature.
	"""
	def __init__(self,term,p_value,N,n,K,genes,pc,es=None):
		mHGTermResult.__init__(self,term,p_value,N,n,K,genes)
		self.pc = pc
		self.enrichment_score = es

	def __repr__(self):
		return "<GOPCASignature: %s / PC %d (p=%.1e; fe=%.1x; es=%.1fx; X=%d; L=%d; %d/%d@%d/%d), gene set %d>" \
				%(self.term.id,self.pc,self.p_value,self.fold_enrichment,self.enrichment_score,self.X,self.L,self.k,self.K,self.n,self.N,hash(self.genes))

	def __str__(self):
		return "<GOPCASignature: %s [%d:%d/%d@%d/%d] (p-value = %.1e, fold enr. = %.1fx, enr. score = %.1fx)>" \
				%(str(self.term),self.pc,self.k,self.K,self.n,self.N,self.p_value,self.fold_enrichment,self.enrichment_socre)

	def __hash__(self):
		return hash(repr(self))

	def __eq__(self,other):
		if type(self) != type(other):
			return False
		elif repr(self) == repr(other):
			return True
		else:
			return False

	def get_pretty_format(self,GO,omit_acc=False,nitty_gritty=True,max_name_length=0):
		term = GO.terms[self.term[0]]
		goterm_genes = GO.get_goterm_genes(term.id)
		details = ''
		if nitty_gritty:
			details = ' [pval=%.1e,fe=%.1fx,E=%.1fx,pc=%d,%d/%d@%d]' \
					%(self.p_value,self.fold_enrichment,self.enrichment_score,self.pc,len(self.genes),self.K,self.n)
			#details = ' [pval=%.1e,fe=%.1fx,pc=%d,%d/%d@%d]' \
			#		%(self.p_value,self.fold_enrichment,self.pc,len(self.genes),self.K,self.n)
		return '%s%s' %(term.get_pretty_format(omit_acc=omit_acc,max_name_length=max_name_length),details)

	@staticmethod
	def from_mHGTermResult(result,pc,es=None):
		return GOPCASignature(result.term,result.p_value,result.N,result.n,result.K,result.genes,pc,es)

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
