cimport cython

import numpy as np
cimport numpy as np

np.import_array()

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

cdef double get_enrichment_score(unsigned char[::1] v, int N, int K, int X, int L, double p_max):
	# calculates enrichment score

	if K == 0 or K == N or K < X:
		return 1.0

	cdef int k = 0
	cdef long double p = 1.0
	cdef int threshold = 0
	cdef float fe
	cdef float e = 1.0
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
			# calculate fold-enrichment only if enough elements have been seen and p-value is small enough
			if k >= X and p <= p_max:
 				fe = (<double>k) / ((<double>K) * ((<double>n) / (<double>N)))
				e = max(e,fe)

	return e
