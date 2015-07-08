#!/usr/bin/env python

import sys
import os
import argparse
import csv
import time

from math import ceil
import itertools as it
import functools as ft
from multiprocessing import Pool

import numpy as np

from sklearn.decomposition import RandomizedPCA

# allow explicit relative imports in executable script
# source: http://stackoverflow.com/a/6655098
if __name__ == '__main__' and __package__ is None:
	parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	sys.path.insert(1, parent_dir)
	import gopca
	__package__ = 'gopca'

from gopca import common
from gopca.tools import misc

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-e','--expression-file',required=True)
	parser.add_argument('-o','--output-file',required=True)

	parser.add_argument('-s','--seed',type=int,default=-1)

	parser.add_argument('-d','--principal-components',type=int,default=20)
	parser.add_argument('-t','--permutations',type=int,default=10000)
	parser.add_argument('-b','--permutation-genes',type=int,default=100)
	parser.add_argument('-j','--jobs',type=int,default=1)

	return parser.parse_args()

def permutation_testing(job_indices,seeds,E,d,b):
	""" Perform PCA for a certain number of permutations. """

	p,n = E.shape
	t = len(job_indices)
	F = np.zeros([t*b,d],dtype=np.float32)

	# turn on reporting if this is the first job
	report = False
	if perm_indices[0] == 0:
		report = True

	# we produce a matrix of T test statistics
	T = np.empty((b*t,d),dtype=np.float32)

	# do the permutations
	total_runs = t*subsets # for calculating progress
	r = 0
	p,n = E.shape
	for j,idx in enumerate(job_indices):

		if report:
			print '\rAnalyzing permuted data: %.1f%% completed...' %(100*(j/float(t))),; sys.stdout.flush()

		np.random.seed(seeds[idx])

		perm_genes = np.random.choice(p,size=b,replace=False) # sample b genes
		E_perm = E.copy()
		for i in perm_genes:
			np.random.shuffle(E_perm[i,:]) # permute row (gene) of expression matrix

		M = RandomizedPCA(n_components=d,random_state=seeds[idx])
		M.fit(E_perm.T)
		W = M.components_.T

		# regression with E that is already centered
		S = E_perm.T.dot(W).T # use this
		x_mean = np.mean(S,axis=1)
		Xc = S - np.tile(x_mean,(n,1)).T
		for pc in range(d):
			s_xx = np.sum(np.power(Xc[pc,:],2.0))
			for i,idx in enumerate(perm_genes):
				s_xy = np.sum(Xc[pc,:]*E_perm[idx,:])
				beta = s_xy / s_xx
				alpha = 0 - beta*x_mean[pc]
				s_sq = np.sum(np.power(E_perm[idx,:]-alpha-beta*S[pc,:],2.0))/float(n-2) # residual variance
				T[j*b+i,pc] = pow(pow(beta,2.0)/(s_sq/s_xx),0.5)

	if report: print 'done!'; sys.stdout.flush()
	return T
	
def store_result(T_sub,b,job_indices,T):
	start = job_indices
	for j,idx in enumerate(job_indices):
		start = idx*b
		stop = start + b
		T[start:stop,:] = T_sub[(j*b):(j*b+b),:]

def main(args):

	expr = args.expression_file
	output_file = args.output_file

	t = args.permutations
	b = args.permutation_genes
	d = args.principal_components
	seed = args.seed
	jobs = args.jobs

	# checks
	assert (t % b) == 0

	# select/generate seed for random number generator
	max_int = np.iinfo(np.int32).max
	if seed < 0:
		seed = np.random.randint(0,max_int)
	np.random.seed(seed)
	print "Seed used:",seed; sys.stdout.flush()

	# read expression
	genes,samples,E = common.read_expression(expr)

	# center expression matrix
	n = E.shape[1]
	E = E - np.tile(np.mean(E,axis=1),(n,1)).T

	# check if there are genes without expression
	sel = np.nonzero(np.amax(E,axis=1)==0)[0]
	if sel.size > 0:
		print "Warning: %d genes with all zero expression values." %(sel.size); sys.stdout.flush()

	# truncate for kicks
	#genes = genes[:500]
	#E = E[:500,:]

	# split up subsamples into distinct jobs
	start = 0
	total_jobs = t/b
	step_size = int(ceil(total_jobs/float(jobs)))
	job_indices = []
	while start < total_jobs:
		job_indices.append(range(start,min(start+step_size,total_jobs),1))
		start += step_size
	#assert len(job_indices) == jobs
	assert sum(len(job_indices[i]) for i in range(jobs)) == total_jobs
	workers = len(job_indices)

	# generate seeds for jobs
	seeds = []
	max_int = np.iinfo(np.int32).max
	for i in range(total_jobs):
		seeds.append(np.random.randint(0,max_int))

	pool = Pool(jobs)

	# perform PCA on unpermuted data
	M = RandomizedPCA(n_components=d)
	M.fit(E.T)
	#W = M.components_.T
	explained = M.explained_variance_ratio_
	print "Explained variance per PC:"
	print explained
	print np.cumsum(explained)
	sys.stdout.flush()

	# run jobs in parallel
	print "Queuing all jobs...", ; sys.stdout.flush()
	p = E.shape[0]
	T_perm = np.zeros((t,d),dtype=np.float32) # for storing output
	pool = Pool(jobs)
	for i in range(jobs):
		func = ft.partial(store_result,b,job_indices=job_indices[i],T=T_perm)
		result = pool.apply_async(permutation_PCA,args=(perm_ind,seeds,E,d,perc),callback=func)
	print "done!"; sys.stdout.flush()

	# waiting for jobs to finish
	t0 = time.time()
	pool.close()
	pool.join()
	t1 = time.time()
	print "All jobs completed (runtime=%.1fs)!" %(t1-t0); sys.stdout.flush()

	# write output
	#labels = ['PC_%d' %(i+1) for i in range(d)]
	#common.write_gene_data(output_file,genes,labels,Z)

	return 0

def run_from_cmdline():
	args = read_args_from_cmdline() 
	return_code = main(args)
	sys.exit(return_code)

if __name__ == '__main__':
	run_from_cmdline()
