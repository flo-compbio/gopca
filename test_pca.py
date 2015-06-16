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

	parser.add_argument('--most-variable-genes',type=int,default=0)
	parser.add_argument('-s','--seed',type=int,default=-1)

	parser.add_argument('-d','--principal-components',type=int,default=20)
	parser.add_argument('-t','--permutations',type=int,default=15)
	parser.add_argument('-p','--permutation-percent',type=float,default=1.0)
	parser.add_argument('-j','--jobs',type=int,default=1)

	return parser.parse_args()

def permutation_PCA(perm_indices,seeds,E,d,perc):
	""" Perform PCA for a certain number of permutations. """

	p,n = E.shape
	t = len(perm_indices)
	C = np.zeros([p,d,t],dtype=np.float64)

	# turn on reporting if this is the first job
	report = False
	if perm_indices[0] == 0:
		report = True

	# define gene index ranges for individual permutations
	start = 0
	p = E.shape[0]
	step_size = int(ceil(p*(perc/100.0)))
	subset_ranges = []
	while start < p:
		stop = min(start + step_size,p)
		subset_ranges.append((start,stop))
		start += step_size

	# do the permutations
	subsets = len(subset_ranges)
	total_runs = t*subsets # for calculating progress
	r = 0
	p,n = E.shape
	for j,idx in enumerate(perm_indices):

		np.random.seed(seeds[idx])
		perm_genes = np.random.permutation(p) # randomly group genes in each permutation cycle

		for k in range(subsets):

			if report:
				print '\rAnalyzing permuted data: %.1f%% completed...' %(100*(r/float(total_runs))),; sys.stdout.flush()

			# figure out which genes to permute
			start,stop = subset_ranges[k]
			E_perm = E.copy()
			sel_genes = perm_genes[start:stop]

			# randomly permute sample labels for subset of genes
			perm_samples = np.random.permutation(n)
			for i in sel_genes:
				E_perm[i,:] = E_perm[i,perm_samples]

			# calculate PCA
			M = RandomizedPCA(n_components=d)
			M.fit(E_perm.T)

			# Important: take the absolute loadings!
			C[sel_genes,:,j] = np.absolute(M.components_[:,sel_genes].T)
			r += 1


	if report: print 'done!'; sys.stdout.flush()
	return C
	
def store_result(C_sub,job_id,perm_indices,C):
	for j,idx in enumerate(perm_indices):
		C[:,:,idx] = C_sub[:,:,j]

def main(args):

	expr = args.expression_file
	output_file = args.output_file

	most_variable_genes = args.most_variable_genes
	t = args.permutations
	d = args.principal_components
	perc = args.permutation_percent
	seed = args.seed
	jobs = args.jobs

	# select/generate seed for random number generator
	max_int = np.iinfo(np.int32).max
	if seed < 0:
		seed = np.random.randint(0,max_int)
	np.random.seed(seed)
	print "Seed used:",seed; sys.stdout.flush()

	# read expression
	genes,samples,E = common.read_expression(expr)

	# select genes with expression
	sel = np.nonzero(np.amax(E,axis=1)>0)[0]
	print "%d genes with nonzero expression." %(sel.size); sys.stdout.flush()
	E = E[sel,:]
	genes = [genes[i] for i in sel]

	# filter genes based on variance
	total_var = np.sum(np.var(E,axis=1,ddof=1))
	p = E.shape[0]
	if most_variable_genes > 0:
		full_var = total_var
		var = np.var(E,axis=1,ddof=1)
		a = np.argsort(var)[::-1]
		sel = np.zeros(p,dtype=np.bool_)
		sel[a[:most_variable_genes]] = True
		sel = np.nonzero(sel)[0]
		E = E[sel,:]
		genes = [genes[i] for i in sel]
		total_var = np.sum(np.var(E,axis=1))
		print "Selected %d most variable genes (removing %.1f%% of the total variance)." \
				%(most_variable_genes,100-100*(total_var/full_var)); sys.stdout.flush()
	
	# truncate for kicks
	#genes = genes[:500]
	#E = E[:500,:]

	# split up subsamples into distinct jobs
	start = 0
	step_size = int(ceil(t/float(jobs)))
	job_indices = []
	while start < t:
		job_indices.append(range(start,min(start+step_size,t),1))
		start += step_size
	assert len(job_indices) == jobs
	assert sum(len(job_indices[i]) for i in range(jobs)) == t

	# generate seeds for jobs
	seeds = []
	max_int = np.iinfo(np.int32).max
	for i in range(t):
		seeds.append(np.random.randint(0,max_int))

	p = E.shape[0]
	C_perm = np.zeros((p,t),dtype=np.float64)
	pool = Pool(jobs)

	# perform PCA on unpermuted data
	M = RandomizedPCA(n_components=d)
	M.fit(E.T)
	C = M.components_.T
	explained = M.explained_variance_ratio_
	print "Explained variance per PC:"
	print explained
	print np.cumsum(explained)
	sys.stdout.flush()

	# run jobs in parallel
	print "Queuing all jobs...", ; sys.stdout.flush()
	p = E.shape[0]
	C_perm = np.zeros((p,d,t),dtype=np.float64) # stores output
	pool = Pool(jobs)
	for i in range(jobs):
		perm_ind = job_indices[i]
		func = ft.partial(store_result,job_id=i,perm_indices=perm_ind,C=C_perm)
		result = pool.apply_async(permutation_PCA,args=(perm_ind,seeds,E,d,perc),callback=func)
	print "done!"; sys.stdout.flush()

	# waiting for jobs to finish
	t0 = time.time()
	pool.close()
	pool.join()
	t1 = time.time()
	print "All jobs completed (runtime=%.1fs)!" %(t1-t0); sys.stdout.flush()

	# calculate Z-scores
	mean = np.mean(C_perm,axis=-1)
	std = np.std(C_perm,axis=-1,ddof=1)
	Z = (np.absolute(C)-mean)/std
	print np.sum(Z>3.0,axis=0)

	# write output
	labels = ['PC_%d' %(i+1) for i in range(d)]
	common.write_gene_data(output_file,genes,labels,Z)

	return 0

if __name__ == '__main__':
	return_code = main(read_args_from_cmdline())
	sys.exit(return_code)
