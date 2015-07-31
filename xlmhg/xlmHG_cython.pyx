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

#cython: profile=False
#cython: wraparound=False
#cython: boundscheck=False
#cython: cdivision=True

cimport cython

import numpy as np
cimport numpy as np

np.import_array()

import sys

cdef extern from "math.h":
	long double fabsl(long double x)
	double NAN

cdef inline int is_equal(long double a, long double b, long double tol):
	# tests equality of two floating point numbers (of type long doube => 80-bit extended precision)
	if a == b:
		return 1
	elif fabsl(a-b)/max(fabsl(a),fabsl(b)) < tol:
		return 1
	else:
		return 0


cdef long double get_hypergeometric_pvalue(\
		long double p, int k, int N, int K, int n):
	# calculates hypergeometric p-value when P(k | N,K,n) is already known
	cdef long double pval = p
	while k < min(K,n):
		p *= (<long double>((n-k)*(K-k)) /\
				<long double>((k+1)*(N-K-n+k+1)))
		pval += p
		k += 1
	return pval


cdef int get_mHG_test_statistic(unsigned char[::1] v, int N, int K, int X, int L,
		long double[::1] mHG_array,
		long double tol):
	# calculates XL-mHG test statistic
	# stores statistic in supplied array, and returns threshold at which minimum is achieved

	if K == 0 or K == N or K < X:
		mHG_array[0] = 1.0
		return 0

	cdef int k = 0
	cdef long double p = 1.0
	cdef long double pval
	cdef long double mHG = 1.1
	cdef int threshold = 0
	cdef int n
	for n in range(L):
		if v[n] == 0:
			# calculate P(k | N,K,n+1) from P(k | N,K,n)
			p *= (<long double>((n+1)*(N-K-n+k)) /\
					<long double>((N-n)*(n-k+1)))
		else:
			# hit one => calculate hypergeometric p-value
			# calculate P(k+1 | N,K,n+1) from P(k | N,K,n)
			p *= (<long double>((n+1)*(K-k)) /\
					<long double>((N-n)*(k+1)))
			k += 1
			if k >= X: # calculate p-value only if enough elements have been seen
				pval = get_hypergeometric_pvalue(p,k,N,K,n+1)

				if pval < mHG and (not is_equal(pval,mHG,tol)):
					# make sure we don't set mHG to something negative
					if pval < 0:
						mHG = 0
					else:
						mHG = pval
					threshold = n+1

	if threshold == 0: # there were not enough positives in v[:L]
		mHG_array[0] = 1.0
	else:
		mHG_array[0] = mHG
	return threshold


cdef long double get_mHG_pvalue(long double s, int N, int K, int X, int L,\
		long double[:,::1] matrix,\
		long double tol):
	# calculates XL-mHG p-value

	# cheap checks
	if s > 1.0 or is_equal(s,1.0,tol):
		return 1.0
	elif s == 0:
		return 0
	elif K == 0 or K >= N or K < X:
		return 0
	elif L > N:
		return 0

	# initialization
	cdef int W = N-K
	cdef int n,k,w
	
	cdef long double p_start = 1.0
	cdef long double p
	cdef long double pval
	matrix[0,0] = 1.0

	# go over all thresholds, except last
	for n in range(1,N):

		if K >= n:
			k = n
			p_start *= ((<long double>(K - n + 1)) /\
					(<long double>(N - n + 1)))
		else:
			k = K
			p_start *= ((<long double>n) /\
					<long double>(n - K))

		if p_start <= 0:
			# not enough floating point precision to calculate p-value
			return <long double>NAN

		p = p_start
		pval = p_start
		w = n - k

		# R is the space of configurations with mHG better than or equal to the one observed
		# - go over all configurations for threshold n
		# - start with highest possible enrichment and then go down
		# - as long as we're in R, all paths going through this configuration are "doomed"
		# - because we're using (K x W) grid instead of parallelogram, "going down" becomes going down and right...

		# no configuration with threshold > L or threshold < X can be in R 
		if n <= L and n >= X:
			# find the first configuration that's not in R
			# this happens when either k < X, or hypergeometric p-value > mHG
			# if k == 0 or w == W, we have hypergeometric p-value = 1
			# since mHG < 1, as soon as k == 0 or w == W, we have left R
			while k >= X and w < W and (pval < s or is_equal(pval,s,tol)):
				# k > 0 is implied
				matrix[k,w] = 0 # we're still in R
				p *= ((<long double>(k*(N-K-n+k))) / (<long double>((n-k+1)*(K-k+1))))
				pval += p
				w += 1
				k -= 1

		# fill in rest of the matrix based on entries for threshold n-1
		while k >= 0 and w <= W:
			if w > 0 and k > 0:
				matrix[k,w] = matrix[k,w-1]*(<long double>(W-w+1))/(<long double>(N-n+1)) +\
					matrix[k-1,w]*(<long double>(K-k+1))/(<long double>(N-n+1))
			elif w > 0:
				matrix[k,w] = matrix[k,w-1]*(<long double>(W-w+1))/(<long double>(N-n+1))

			elif k > 0:
				matrix[k,w] = matrix[k-1,w]*(<long double>(K-k+1))/(<long double>(N-n+1))

			w += 1
			k -= 1

	return 1.0 - (matrix[K,W-1] + matrix[K-1,W])

def xlmhg(unsigned char[::1] v, int N, int K, int X, int L, int use_upper_bound, long double[:,::1] matrix, double tol, double pval_thresh):

	cdef int n
	cdef long double[::1] s_array = np.zeros(1,dtype=np.longdouble)
	cdef long double s, pval
	cdef double s_double, pval_double
	cdef long double tol_ld = <long double>tol
	cdef long double pval_thresh_ld = <long double>pval_thresh

	# get XL-mHG test statistic and corresponding threshold
	n = get_mHG_test_statistic(v, N, K, X, L, s_array, tol_ld)
	s = s_array[0]
	if s > 1.0 or is_equal(s,1.0,tol_ld):
		return n,1.0,1.0

	# get XL-mHG p-value (either exact or using upper bound)
	pval_double = 0
	if ((not is_equal(s,pval_thresh_ld,tol_ld)) and s > pval_thresh_ld) or use_upper_bound != 0:
		# use upper bound
		pval_double = <double>(min(1.0,s*<long double>K))
	else:
		# calculate exact p-value
		pval = get_mHG_pvalue(s, N, K, X, L, matrix, tol_ld)
		pval_double = <double>pval

	s_double = <double>s
	return n,s_double,pval_double
