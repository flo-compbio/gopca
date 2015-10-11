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
from sklearn.decomposition import PCA, RandomizedPCA

from gopca import common
from gopca.printf import printf

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-e','--expression-file',required=True)
	parser.add_argument('-D','--components',type=int,required=True)
	parser.add_argument('-o','--output-file',required=True)

	parser.add_argument('-G','--select-variable-genes',type=int,default=0)
	parser.add_argument('-t','--permutations',type=int,default=15)
	parser.add_argument('-z','--zscore-thresh',type=float,default=2.0)

	parser.add_argument('-fs','--figure-size',type=float,help='in inches',nargs=2,default=[18,10])
	parser.add_argument('-fr','--figure-resolution',type=int,help='in dpi',default=150)
	parser.add_argument('-ff','--figure-font-size',type=int,help='in pt',default=24)
	parser.add_argument('-fm','--figure-font-family',default='serif')
	parser.add_argument('-ft','--figure-use-tex',action='store_true')
	parser.add_argument('-b','--matplotlib-backend',default=None)

	parser.add_argument('-s','--seed',type=int,default=None)
	parser.add_argument('--quiet',action='store_true')

	return parser.parse_args()

def message(m,quiet,flush=True,endline=True):
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
	n_components = args.components
	output_file = args.output_file

	sel_var_genes = args.select_variable_genes
	t = args.permutations
	zscore_thresh = args.zscore_thresh
	seed = args.seed
	quiet = args.quiet

	# figure size
	fig_size = args.figure_size
	fig_res = args.figure_resolution

	# figure text
	fig_use_tex = args.figure_use_tex
	fig_font_size = args.figure_font_size
	fig_font_family = args.figure_font_family

	mpl_backend = args.matplotlib_backend
	# set seed for random number generator
	if seed is None:
		seed = np.random.randint(int(1e9))
	np.random.seed(seed)
	message('Using seed: %d' %(seed),quiet)

	# checks
	assert os.path.isfile(expression_file)
	assert t >= 2
	assert zscore_thresh >= 0

	# read expression
	genes,samples,E = common.read_expression(expression_file)

	# filter for most variable genes
	if sel_var_genes > 0:
		var = np.var(E,axis=1)
		a = np.argsort(var)
		a = a[::-1]
		genes = [genes[i] for i in a[:sel_var_genes]]
		E = E[a[:sel_var_genes]]
		
	message('Expression matrix shape: ' + str(E.shape),quiet)
	#E += (np.random.rand(*E.shape)*1e-4)

	# do PCA on unpermuted data
	message('Performing PCA on unpermuted data...',quiet=quiet,endline=False)
	p,n = E.shape
	n_comps = min(min(p,n-1),n_components)
	M = PCA(n_components = n_comps)
	M.fit(E.T)
	d = M.explained_variance_ratio_.copy()
	message('done!',quiet=quiet)

	# do permutations
	message('Performing permutations...',quiet=quiet,endline=False)
	thresh = common.get_pc_explained_variance_threshold(E,zscore_thresh,t,seed=seed)
	message('done!',quiet=quiet)

	d_est = np.sum(d >= thresh)
	args.result = d_est
	message('Number of principal components with z-score >= %.1f: %d' \
			%(zscore_thresh,d_est),quiet)

	# plotting
	message('Plotting...',quiet=quiet,endline=False)

	import matplotlib as mpl
	if mpl_backend is not None:
		mpl.use(mpl_backend)
	import matplotlib.pyplot as plt
	from matplotlib import rc

	if fig_use_tex: rc('text',usetex=True)
	rc('font',family=fig_font_family,size=fig_font_size)
	rc('figure',figsize=(fig_size[0],fig_size[1]))
	rc('savefig',dpi=fig_res)

	b1 = None
	b2 = None
	if d_est > 0:
		b1 = plt.bar(np.arange(d_est)-0.4,100*d[:d_est],width=0.8,edgecolor='none',color='gold')
	if n_comps > d_est:
		b2 = plt.bar(np.arange(d_est,n_comps)-0.4,100*d[d_est:],width=0.8,edgecolor='none',color='skyblue')
	#plt.plot([-0.5,n_comps-0.5],[mean_null,mean_null],color='pink',lw=1.0)
	#plt.plot([-0.5,n_comps-0.5],[mean_null+std_null,mean_null+std_null],'--',color='red',lw=1.0)
	plt.plot([-0.5,n_comps-0.5],[100*thresh,100*thresh],':',color='gray',lw=2.0)

	if fig_use_tex:
		plt.ylabel(r'Fraction of Variance Explained (\%)')
	else:
		plt.ylabel('Fraction of Variance Explained (%)')

	common.simpleaxis(plt.gca())

	plt.xticks(np.arange(n_comps),np.arange(n_comps)+1)
	plt.xlabel('Principal Component')
	plt.xlim(-0.5,n_comps-0.5)

	ax2 = plt.twinx()
	plt.sca(ax2)
	plt.plot(np.arange(n_comps),100*np.cumsum(d),color='purple',lw=2.0,marker='o',markerfacecolor='purple',markeredgecolor='none',markersize=10)

	plt.yticks(color='purple')
	if fig_use_tex:
		plt.ylabel(r'Cumulative Fraction of Variance Explained (\%)',color='purple')
	else:
		plt.ylabel('Cumulative Fraction of Variance Explained (%)',color='purple')

	plt.gca().spines['top'].set_visible(False)
	#common.simpleaxis(plt.gca())


	message('Saving to file...',quiet,endline=False)
	#plt.gcf().set_size_inches(fig_dim[0],fig_dim[1])
	plt.savefig(output_file,bbox_inches='tight')
	message('done!',quiet)

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
