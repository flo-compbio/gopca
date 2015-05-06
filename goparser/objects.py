# Copyright (c) 2015 Florian Wagner
#
# This file is part of GOParser.
#
# GOParser is free software: you can redistribute it and/or modify
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

class GOTerm(object):

	short_ns = {'biological_process': 'BP', 'molecular_function': 'MF', 'cellular_component': 'CC'}

	def __init__(self,id_,name,namespace,is_a,part_of,definition=None):
		self.id = id_
		self.name = name
		self.namespace = namespace
		self.definition = definition

		# to store immediate parents/wholes
		self.is_a = is_a.copy()
		self.part_of = part_of.copy()

		# to store immediate children/parts
		self.children = set()
		self.parts = set()

		# to store all descendants/ancestors
		self.descendants = None
		self.ancestors = None

	def __repr__(self):
		return "<GOTerm %s>" %(self.id)
	def __str__(self):
		return "<GOTerm: %s>" %(self.get_pretty_format())

	def __eq__(self,other):
		if type(self) != type(other):
			return False

		elif self is other:
			return True
			
		elif self.id == other.id:
			return True

		else:
			return False

	def __hash__(self):
		return hash(repr(self))

	def get_namespace_short(self):
		return self.short_ns[self.namespace]

	def get_id(self):
		return self.id

	def get_acc(self): # gets the accession number as integer
		return int(self.id[3:])

	def get_pretty_format(self,omit_acc=False,max_name_length=0):
		#print self.namespace
		name = self.name
		if max_name_length >= 3 and len(name) > max_name_length:
			name = name[:(max_name_length-3)] + '...'
		if omit_acc: return "%s: %s" %(self.short_ns[self.namespace], name)
		#else: return "%s: %s (GO:%07d)" %(self.short_ns[self.namespace], self.name, self.acc)
		else: return "%s: %s (%s)" %(self.short_ns[self.namespace], name, self.id)


class GOAnnotation(object):

	#uniprot_pattern = re.compile("([A-Z][A-Z0-9]{5})(?:-(\d+))?")
	#short_ns = {'biological_process': 'BP', 'molecular_function': 'MF', 'cellular_component': 'CC'}

	def __init__(self,target,term,evidence,pubmed_id=None,uniprot=None):
		assert target is not None and target != ''
		assert evidence is not None and evidence != ''
		assert type(term) == GOTerm
		self.target = target # a gene name
		self.evidence = evidence # GO evidence code
		self.term = term # a GOTerm object
		self.pubmed_id = pubmed_id # PubMed ID, optional
		self.uniprot = uniprot # Uniprot identifier, optional

	evidence_name = {\
			'EXP': 'experiment',\
			'IDA': 'direct assay',\
			'IPI': 'physical interaction',\
			'IMP': 'mutant phenotype',\
			'IGI': 'genetic interaction',\
			'IEP': 'expression pattern',\
			'ISS': 'sequence or structural similarity',\
			'ISO': 'sequence orthology',\
			'ISA': 'sequence alignment',\
			'ISM': 'sequence model',\
			'IGC': 'genomic context',\
			'IBA': 'biological aspect of ancestor',\
			'IBD': 'biological aspect of descendant',\
			'RCA': 'reviewed computational analysis',\
			'TAS': 'traceable author statement',\
			'NAS': 'non-traceable author statement',\
			'IC' : 'inferred by curator',\
			'ND' : 'no biological data available',\
			'IEA': 'inferred from electronic annotation'\
			}

	evidence_type = {\
			'EXP': 'experimental',\
			'IDA': 'experimental',\
			'IPI': 'experimental',\
			'IMP': 'experimental',\
			'IGI': 'experimental',\
			'IEP': 'experimental',\
			'ISS': 'computational',\
			'ISO': 'computational',\
			'ISA': 'computational',\
			'ISM': 'computational',\
			'IGC': 'computational',\
			'IBA': 'computational',\
			'IBD': 'computational',\
			'RCA': 'computational',\
			'TAS': 'literature',\
			'NAS': 'literature',\
			'IC' : 'curator',\
			'ND' : 'curator',\
			'IEA': 'automatic'\
			}

	evidence_type_short = {\
			'EXP': 'exp.',\
			'IDA': 'exp.',\
			'IPI': 'exp.',\
			'IMP': 'exp.',\
			'IGI': 'exp.',\
			'IEP': 'exp.',\
			'ISS': 'comp.',\
			'ISO': 'comp.',\
			'ISA': 'comp.',\
			'ISM': 'comp.',\
			'IGC': 'comp.',\
			'IBA': 'comp.',\
			'IBD': 'comp.',\
			'RCA': 'comp.',\
			'TAS': 'lit.',\
			'NAS': 'lit.',\
			'IC' : 'cur.',\
			'ND' : 'cur.',\
			'IEA': 'autom.'\
			}

	def __repr__(self):
		uniprot = ''
		if self.uniprot is not None:
			uniprot = self.uniprot
		pmid = ''
		if self.pubmed_id is not None:
			pmid = self.pubmed_id
		return "<GOAnnotation %s>" %(', '.join([self.target,uniprot,repr(self.term),self.evidence,pmid]))

	def __str__(self):
		return "<GOAnnotation: %s>" %(self.get_pretty_format())

	def __eq__(self,other):
		if type(self) != type(other):
			return False

		elif self is other:
			return True

		elif self.target == other.target \
				and self.uniprot == other.uniprot \
				and self.term == other.term \
				and self.evidence == other.evidence \
				and self.pubmed_id == other.pubmed_id:
			return True

		else:
			return False

	def __hash__(self):
		return hash(repr(self))

	def get_pretty_format(self,uniprot_names = None):
		#return '%s\t%s\t%s/%s\t%s' %(short_ns[self.term.namespace], self.term.name, self.evidence_type[self.evidence], self.evidence_name[self.evidence],self.term.acc)
		pmid = ''
		if self.pubmed_id is not None: pmid = self.pubmed_id
		uniprot = ''
		if self.uniprot is not None:
			uniprot = self.uniprot
			if uniprot_names is not None:
				try:
					n = uniprot_names[self.uniprot]
					uniprot = n[0]+'/'+n[1]
				except KeyError:
					pass
		#return '%s\t%s\t%s\t%s/%s\t%s\t%s' %(self.term.get_namespace_short(), self.term.get_id(), self.term.name, self.evidence_type[self.evidence], self.evidence_name[self.evidence],pmid,name)
		pretty = '\t'.join([self.target,uniprot,self.evidence,pmid,self.term.get_pretty_format()])
		return pretty


