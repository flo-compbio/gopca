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



import csv
import numpy as np

from tools import misc

def read_meta(fn):
	meta = {}
	with open(fn) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		for l in reader:
			meta[l[0]] = l[1:]
	return meta

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

def write_expression(ofn,genes,samples,E):
	n,p = E.shape
	with open(ofn,'w') as ofh:
		writer = csv.writer(ofh,dialect='excel-tab',lineterminator='\n',quoting=csv.QUOTE_NONE)
		writer.writerow([''] + samples)
		for i in range(n):
			writer.writerow([genes[i]] + ['%.5f' %(E[i,j]) for j in range(p)])

def get_signature(genes,E,sig_genes):
    m = len(sig_genes)
    n,p = E.shape
    S = np.zeros((m,p),dtype=np.float64)
    for i,g in enumerate(sig_genes):
        idx = misc.bisect_index(genes,g)
        S[i,:] = E[idx,:]
        S[i,:] -= np.mean(S[i,:])
        S[i,:] /= np.std(S[i,:],ddof=1)
    sig = np.mean(S,axis=0)
    return sig

def get_signature_robust(genes,E,sig_genes):
	m = len(sig_genes)
	n,p = E.shape
	S = np.zeros((m,p),dtype=np.float64)
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
