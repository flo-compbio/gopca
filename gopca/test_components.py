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

import numpy as np
from scipy.stats import norm
from sklearn.decomposition import PCA, RandomizedPCA
#from sklearn.decomposition import RandomizedPCA

from gopca import common
from gopca.printf import printf

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='GO-PCA')

	parser.add_argument('-e','--expression-file',required=True)
	parser.add_argument('-t','--permutations',type=int,default=15)
	parser.add_argument('-c','--test-components',type=int,default=50)
	parser.add_argument('-q','--qval-thresh',type=float,default=0.05)

	parser.add_argument('-s','--seed',type=int,default=123456789)
	parser.add_argument('--quiet',action='store_true')

	return parser.parse_args()

def message(m,quiet,flush=True,endline=False):
	if not quiet:
		end = ' '
		if endline:
			end = '\n'
		printf(m,end=end)
		if flush:
			sys.stdout.flush()

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	expression_file = args.expression_file
	t = args.permutations
	test_components = args.test_components
	qval_thresh = args.qval_thresh
	seed = args.seed
	quiet = args.quiet

	# set seed for random number generator
	if seed is not None:
		np.random.seed(seed)

	# checks
	assert os.path.isfile(expression_file)
	assert t >= 2
	assert 0 < qval_thresh <= 1.0

	# read expression
	genes,samples,E = common.read_expression(expression_file)

	#E += (np.random.rand(*E.shape)*1e-4)

	# do PCA on unpermuted data
	p,n = E.shape
	n_comps = test_components
	all_comps = min(p,n-1)
	if all_comps < test_components:
		n_comps = all_comps
		message('Warning: The number of PCs to test that was specified is %d. However, there are only %d principal components!' \
				%(test_components,n_comps),quiet=quiet)
		print
		
	M = PCA(n_components = n_comps)
	M.fit(E.T)
	d = M.explained_variance_ratio_.copy()

	# do permutations
	message('Performing permutations...',quiet=quiet,endline=False)
	d_max_null = np.empty(t,dtype=np.float64)
	E_perm = np.empty((p,n),dtype=np.float64)
	M_null = RandomizedPCA(n_components = 1, random_state=seed)
	for j in xrange(t):
		message('%d...' %(j+1),quiet=quiet,endline=False)

		for i in range(p):
			E_perm[i,:] = E[i,np.random.permutation(n)]

		M_null.fit(E_perm.T)
		d_max_null[j] = M_null.explained_variance_ratio_[0]
	message('done!',quiet=quiet)

	mean_null = np.mean(d_max_null)
	std_null = np.std(d_max_null,ddof=1)

	# calculate z-scores and corresponding p-values
	zscores = (d-mean_null) / std_null
	pvals = norm.sf(zscores)

	if test_components > n_comps:
		pvals = np.r_[pvals,[1.0]*(test_components-n_comps)]
	assert pvals.size == test_components

	# conservatively force monotonicity of p-values
	for i in xrange(1,n_comps):
		pvals[i] = max(pvals[i-1],pvals[i])
	assert np.all(np.argsort(pvals,kind='mergesort') == np.arange(test_components)) # stable sort

	# calculate q-values
	qvals = common.get_qvalues(pvals)
	assert np.all(np.argsort(qvals,kind='mergesort') == np.arange(test_components)) # stable sort

	significant = np.sum(qvals <= qval_thresh)
	args.output = significant
	message('Number of significant principal components at FDR threshold of q=%.3f: %d (p-value cutoff = %.1e)' \
			%(qval_thresh,significant,pvals[significant-1]),quiet)

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
