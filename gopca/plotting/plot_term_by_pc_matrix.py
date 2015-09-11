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
import csv
import cPickle as pickle
from math import ceil
from collections import OrderedDict,Counter

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy import stats 
from sklearn.decomposition import PCA

import matplotlib as mpl
#from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc

# allow explicit relative imports in executable script
# source: http://stackoverflow.com/a/6655098
if __name__ == '__main__' and __package__ is None:
	parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	sys.path.insert(1, parent_dir)
	__package__ = 'gopca'
	#from .. import gopca
	import gopca

from gopca import common
from gopca.tools import misc
from gopca.go_pca_objects import GOPCASignature
#from gopca.xlmhg.xlmHG_cython import mHG_test
import xlmhg

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-g','--gopca-file',required=True)
	parser.add_argument('-a','--annotation-fie',required=True)
	parser.add_argument('-o','--output-file',required=True)

	# visualization options
	parser.add_argument('--pvalue-show-below',type=float,default=1e-4)
	parser.add_argument('--pvalue-best',type=float,default=1e-10)
	parser.add_argument('--pvalue-worst',type=float,default=1.0)
	parser.add_argument('--pvalue-ticks',type=float,nargs='+',default=[4,6,8,10])
	parser.add_argument('-r','--reverse-order',required=False,action='store_true')

	# cosmetics
	parser.add_argument('--term-max-name-length',type=int,default=50)
	parser.add_argument('--term-show-id',action='store_true')

	parser.add_argument('-s','--figure-size',type=int,help='in inches',nargs=2,default=[10,15])
	parser.add_argument('-f','--figure-font-size',type=int,help='in pt',default=24)
	parser.add_argument('-m','--figure-font-family',default='serif')
	parser.add_argument('-r','--figure-resolution',type=float,help='in dpi',default=150.0)
	parser.add_argument('-c','--figure-colormap',default='BuPu')

	parser.add_argument('-t','--use-tex',action='store_true')
	parser.add_argument('-b','--backend',default=None)

	parser.add_argument('--dotcolor',default='yellow')
	parser.add_argument('--dotsize',type=float,default=50)

	return parser.parse_args()

def main(args=None):

	if args is None:
		# read command line arguments
		args = read_args_from_cmdline()

	gopca_file = args.gopca_file
	annotation_file = args.annotation_file
	output_file = args.output_file

	pval_show_below = args.pvalue_show_below
	pval_best = args.pvalue_best
	pval_worst = args.pvalue_worst
	pval_ticks = args.pvalue_ticks
	reverse_order = args.reverse_order

	fig_size = args.figure_size
	fig_font_size = args.figure_font_size
	fig_font_family = args.figure_font_family
	fig_dpi = args.figure_resolution
	fig_cmap = args.figure_colormap
	fig_dotcolor = args.dotcolor
	fi_dotsize = args.dotsize

	max_name_length = args.term_max_name_length
	omit_acc = not args.term_show_id

	# read GO-PCA result
	result = common.read_gopca_result(gopca_file)

	# read annotations
	annotations = common.read_annotations(annotation_file)

	genes = result.genes
	signatures = result.signatures
	W = result.W # loading matrix
	S = result.S # signature matrix
	mHG_X_frac = result.mHG_X_frac
	mHG_X_min = result.mHG_X_min
	mHG_L = result.mHG_L

	# determine max K
	term_genes = dict([sig.term,annotations[sig.term] for sig in signatures])
	K_max = max(len(tg) for tg in term_genes.itervalues())

	# generate signature labels
	labels = [sig.get_label(max_name_length) for sig in signatures]

	# order signatures using hierarchical clustering
	order = common.cluster_signatures(S)
	if reverse_order:
		order = order[::-1]
	S = S[order,:]
	labels = [labels[i] for i in order]
	signatures = [signatures[i] for i in order]

	# test association
	q = len(signatures)
	A = np.zeros((q,2*test_pc),dtype=np.float64)
	p = E.shape[0]
	matrix = np.zeros((K_max+1,p+1),dtype=np.longdouble)

	all_genes = set(genes)
	sig_term_genes = [sig.term[0]

	for pc in range(test_pc):
		a_pos = np.argsort(-W[:,pc])
		a_neg = a_pos[::-1]
		for i,(sig,tg) in enumerate(zip(signatures,term_genes)):
			K = len(tg)
			X = max(mHG_X_min,int(ceil(mHG_X_frac*K)))
		
			v = np.zeros(p,dtype=np.uint8)
			for g in term_genes:
				idx = misc.bisect_index(genes,g)
				v[idx] = 1

			v_sorted = np.ascontiguousarray(v[a_pos])
			threshold,_,pval = mHG_test(v_sorted,p,K,mHG_L,X,mat=matrix)
			A[i,pc*2] = -np.log10(pval)

			v_sorted = np.ascontiguousarray(v[a_neg])
			threshold,_,pval = mHG_test(v_sorted,p,K,mHG_L,X,mat=matrix)
			A[i,pc*2+1] = -np.log10(pval)

	# prepare for plotting
	rc('font',family=fig_font_family,size=fig_font_size)
	rc('figure',figsize=fig_dim)
	rc('savefig',dpi=fig_dpi)

	if use_tex:
		rc('text',usetex=True)
		preamble = mpl.rcParams['text.latex.preamble']
		add = r'\usepackage{bm}'
		if add not in preamble:
			mpl.rcParams['text.latex.preamble'].append(add)

	# for each signature, mark PC that was originally used to generate it
	# plot this first, otherwise colormap gets messed up
	q = len(signatures)
	for i in range(q):
		#sel = np.nonzero(np.absolute(A[i,:])>=-np.log10(go_pvalue_threshold))[0]
		j = (abs(signatures[i].pc)-1)*2
		if signatures[i].pc < 0: j+=1
		#plt.scatter([j],[i],marker='o',facecolor=dotcolor,color='none',zorder=100,s=dotsize)
		plt.scatter([j],[i],color=dotcolor,zorder=100,s=dotsize,marker='x')

	# plot heatmap
	A[np.absolute(A)<-np.log10(pval_show_below)] = np.nan # hide insignificant associations
	vmin = -np.log10(pval_worst)
	vmax = -np.log10(pval_best)
	plt.imshow(A,interpolation='none',vmin=vmin,vmax=vmax,cmap=cmap,zorder=20)

	# plot separation lines
	for pc in range(test_pc-1):
		plt.plot([pc*2+1.5,pc*2+1.5],[-0.5,q-0.5],color='gray',linewidth=1.0,zorder=50)
	plt.plot([-0.5,test_pc*2+1.5],[-0.5,-0.5],color='black',linewidth=2.0,zorder=50) # fix some z-order issues

	# configure axes
	plt.gca().yaxis.set_ticks_position('left')
	plt.gca().xaxis.set_ticks_position('top')
	plt.gca().xaxis.set_label_position('top')
	plt.gca().spines['right'].set_visible(False)
	plt.gca().spines['bottom'].set_visible(False)
 
	#plt.gca().xaxis.tick_top()
	plt.xticks(np.arange(0,test_pc*2,2)+0.5,np.arange(test_pc)+1,size='small')
	plt.xlabel(r'Principal Components',labelpad=10,size='small')
	plt.xlim(-0.5,test_pc*2-0.5)
	plt.gca().tick_params(top='off')

	plt.yticks(np.arange(q),labels,size='x-small')
	plt.ylabel(r'GO Terms',size='small')
	plt.ylim(q-0.5,-0.5)
	#plt.ylim(-0.5,q-0.5)
	plt.grid(which='both',axis='y',zorder=-20) # z-order is ignored here

	# plot colorbar
	cbar = plt.colorbar(pad=0.02,shrink=0.5)
	cbar.set_ticks(pval_ticks)
	#cbar.ax.tick_params(labelsize='medium')
	cbar.set_label(r"$\bm{-\log_{10}}$ p-value")

	plt.show()

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
