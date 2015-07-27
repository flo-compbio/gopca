import numpy as np

def get_enrichment_score(v, N, K, X, L, p_max):
	# calculates enrichment score

	if K == 0 or K == N or K < X:
		return 1.0

	k = 0
	p = 1.0
	e = 1.0
	for n in range(L):
		if v[n] == 0:
			# calculate P(k | N,K,n+1) from P(k | N,K,n)
			p *= (float((n+1)*(N-K-n+k)) /\
					float((N-n)*(n-k+1)))
		else:
			# hit one => calculate hypergeometric p-value
			# calculate P(k+1 | N,K,n+1) from P(k | N,K,n)
			p *= (float((n+1)*(K-k)) /\
					float((N-n)*(k+1)))
			k += 1
			# calculate fold-enrichment only if enough elements have been seen and p-value is small enough
			if k >= X and p <= p_max:
 				fe = k / (K * (n / float(N)))
				e = max(e,fe)

	return e
