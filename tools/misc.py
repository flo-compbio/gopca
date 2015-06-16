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
import bisect

def flatten(l):
        return [item for sublist in l for item in sublist] # incomprehensible list comprehension

def bisect_index(a, x):
	'Locate the leftmost value exactly equal to x'
	i = bisect.bisect_left(a, x)
	if i != len(a) and a[i] == x:
		return i
	raise ValueError

def argsort(seq):
	#http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
	#by unutbu
	return sorted(range(len(seq)), key=seq.__getitem__)

def argmin(seq):
	return argsort(seq)[0]

def argmax(seq):
	return argsort(seq)[-1]

def read_single(fn):
	data = []
	with open(fn) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		for l in reader:
			data.append(l[0])
	return data

def read_all(fn,m='r'):
	data = []
	with open(fn,m) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		for l in reader:
			data.append(l)
	return data

def read_chromosome_lengths(fn):
	chromlen = {}
	with open(fn) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		for l in reader:
			chromlen[l[0]] = int(l[1])
	return chromlen

def read_chromlen(fn):
	return read_chromosome_lengths(fn)

def read_enrichment(fn):
	return read_all_columns(fn)

def read_genes(fn):
	return read_single(fn)

def read_genes_and_chromosomes(fn):
	genes = []
	chromosomes = []
	with open(fn) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		i = 0
		for l in reader:
			genes.append(l[0])
			chromosomes.append(set(l[1].split(',')))
			i += 1
	return genes,chromosomes

def read_goterms(fn):
	return read_single_column(fn)

def read_terms(fn):
	return read_all_columns(fn)


