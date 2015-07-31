# Copyright (c) 2015 Florian Wagner
#
# This file is part of the Python/Cython implementation of the XL-mHG.
#
# The Python/Cython implementation of the XL-mHG
# is free software: you can redistribute it and/or modify
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

from math import isnan

import numpy as np

from xlmHG_cython import xlmhg

def mHG_test(v, X, L, K=None, mat=None, use_upper_bound=False, verbose=False, tol=1e-16, pval_thresh=1.0, perform_checks=True):
	"""
	Front-end for the XL-mHG test.
	"""

	### type checks
	if perform_checks:
		assert type(use_upper_bound) == bool or type(use_upper_bound) == int

		# check vector
		assert type(v) == np.ndarray
		assert v.ndim == 1
		assert v.dtype == np.uint8
		if not v.flags.c_contiguous:
			print >> sys.stderr, 'Warning: mHG_test called with vector ("v" parameter) that is not C-contiguous!'
			v = np.ascontiguousarray(v)

		# check parameters
		assert type(X) == int
		assert type(L) == int
		if K is not None:
			assert type(K) == int
		assert type(tol) == float
		assert type(pval_thresh) == float

		# check matrix
		if mat is not None:
			assert type(mat) == np.ndarray
			assert mat.dtype == np.longdouble
			assert mat.ndim == 2
			if not mat.flags.c_contiguous:
				print >> sys.stderr, 'Warning: mHG_test called with matrix ("mat" parameter) that is not C-contiguous!'
				mat = np.ascontiguousarray(mat)

	N = v.size
	# determine K (if not supplied)
	if K is None:
		K = np.nonzero(v)[0].size

	# allocate matrix (if not supplied)
	if mat is None:
		mat = np.empty((K+1,N-K+1),dtype=np.longdouble)

	use_upper_bound == int(use_upper_bound)

	### sanity checks
	if perform_checks:
		assert mat.shape[0] >= K+1 and mat.shape[1] >= N-K+1
		assert N >= 0
		assert 0 <= K <= N
		assert 0 <= L <= N
		assert 0 <= X <= K

	# special cases
	if K == 0 or K == N: # check if we have any positives at all, or if all entries are positives
		return 0,1.0,1.0

	n,s,pval = xlmhg(v,N,K,X,L,use_upper_bound,mat,tol,pval_thresh)

	# check whether floating point accuracy was insufficient for calculation of the p-value
	if isnan(pval) or pval <= 0:
		# if so, use upper bound instead
		pval = min(1.0,s*K)

	return n,s,pval
