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


import os
import csv
import numpy as np

from genometools import misc

def read_meta(fn):
	meta = {}
	with open(fn) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		for l in reader:
			meta[l[0]] = l[1:]
	return meta

def read_gene_data(fn,case_insensitive=False):
	labels = None
	genes = []
	data = []
	with open(fn) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		labels = reader.next()[1:]
		for l in reader:
			g = l[0]
			if case_insensitive: g = g.upper()
			genes.append(g)
			data.append(l[1:])
	D = np.float64(data)
	return genes,labels,D

def read_expression(fn,case_insensitive=False):
	samples = None
	genes = []
	expr = []
	with open(fn) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		samples = reader.next()[1:]
		for l in reader:
			g = l[0]
			if case_insensitive: g = g.upper()
			genes.append(g)
			expr.append(l[1:])
	E = np.float64(expr)
	return genes,samples,E

def write_expression(output_file,genes,samples,E):
	write_gene_data(output_file,genes,samples,E)

def write_gene_data(output_file,genes,labels,D):
	p = len(genes)
	n = len(labels)
	assert D.shape == (p,n)
	with open(output_file,'w') as ofh:
		writer = csv.writer(ofh,dialect='excel-tab',lineterminator=os.linesep,quoting=csv.QUOTE_NONE)
		writer.writerow(['.'] + labels)
		for i,g in enumerate(genes):
			writer.writerow([g] + ['%.5f' %(D[i,j]) for j in range(n)])

def get_signature_expression(genes,E,sig_genes):
	p_sig = len(sig_genes)
	p,n = E.shape
	S = np.zeros((p_sig,n),dtype=np.float64)
	for i,g in enumerate(sig_genes):
		idx = misc.bisect_index(genes,g)
		S[i,:] = E[idx,:]
		S[i,:] -= np.mean(S[i,:])
		S[i,:] /= np.std(S[i,:],ddof=1)
	sig = np.mean(S,axis=0)
	return sig

def get_signature_expression_robust(genes,E,sig_genes):
	p_sig = len(sig_genes)
	p,n = E.shape
	S = np.zeros((p_sig,n),dtype=np.float64)
	for i,g in enumerate(sig_genes):
		idx = misc.bisect_index(genes,g)
		S[i,:] = E[idx,:]
		med = np.median(S[i,:])
		mad = np.median(np.absolute(S[i,:]-med))
		std = 1.4826*mad
		S[i,:] -= med
		S[i,:] /= std
	sig = np.mean(S,axis=0)
	return sig

def get_signature_label(GO,sig,max_length=40):
	count = ' (%d:%d/%d)' %(sig.pc,len(sig.genes),sig.K)
	enr = sig.enrichment
	return GO.terms[enr.term[0]].get_pretty_format(omit_acc=True,max_name_length=max_length) + count

def variance_filter(genes,E,top):
	# filter genes by variance
	a = np.argsort(np.var(E,axis=1,ddof=1))[::-1]
	n = E.shape[0]
	sel = np.zeros(n,dtype=np.bool_)
	sel[a[:top]] = True
	sel = np.nonzero(sel)[0]
	genes = [genes[i] for i in sel]
	E = E[sel,:]
	return genes,E
