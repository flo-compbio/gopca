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

class mHGTermResultWithPC(mHGTermResult):
	"""
	Stores mHG result for one particular term.
	"""
	def __init__(self,term,p_value,N,n,K,genes,pc):
		mHGTermResult.__init__(self,term,p_value,N,n,K,genes)
		self.pc = pc

	def __repr__(self):
		return "<mHGTermEnrichment of term '%s', PC %d, %d genes (hash:%d)>" %(self.term[0],self.pc,len(self.genes),hash(self.genes))

	def __str__(self):
		return "<mHG_Enrichment of term '%s' in PC %d:  p-value = %.1e, fold enrichment = %.2fx, %d genes>" \
				%(str(self.term),self.pc,self.p_value,self.fold_enrichment, len(self.genes))

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
			details = ' [p=%.1e,e=%.2fx,c=%d,%d/%d@%d]' %(self.p_value,self.fold_enrichment,self.pc,len(self.genes),self.K,self.n)
		return '%s%s' %(term.get_pretty_format(omit_acc=omit_acc,max_name_length=max_name_length),details)

	@staticmethod
	def from_mHGTermResult(pc,result):
		return mHGTermResultWithPC(result.term,result.p_value,result.N,result.n,result.K,result.genes,pc)

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
