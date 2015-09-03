#!/usr/bin/env python2.7

import sys
import os
import argparse
import csv
import cPickle as pickle
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
from gopca import fdr
from gopca.tools import misc
from gopca.go_pca_objects import mHGTermResultWithPC
from gopca.xlmhg.xlmHG_cython import mHG_test

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-e','--expression-file',required=True)
	parser.add_argument('-r','--result-file',required=True)
	parser.add_argument('-g','--go-pickle-file',required=True)

	# GO-PCA parameters
	#parser.add_argument('-X','--go-mHG-X',type=int,required=True)
	#parser.add_argument('-L','--go-mHG-L',type=int,required=True)
	#parser.add_argument('-p','--go-pvalue-threshold',type=float,default=1e-6)
	#parser.add_argument('-q','--go-enrichment-fdr',type=float,default=0.05) # ignored if no z-score file is provided (-z)

	# visualization options
	parser.add_argument('--pvalue-show-below',type=float,default=1e-4)
	parser.add_argument('--pvalue-best',type=float,default=1e-10)
	parser.add_argument('--pvalue-worst',type=float,default=1.0)
	parser.add_argument('--pvalue-ticks',type=float,nargs='+',default=[4,6,8,10])
	parser.add_argument('-r','--reverse-order',required=False,action='store_true')

	#parser.add_argument('-o','--output-file',required=True)

	# cosmetics
	parser.add_argument('--GO-term-max-name-length',type=int,default=50)
	parser.add_argument('--GO-term-show-id',action='store_true')

	parser.add_argument('-d','--figure-dimensions',type=int,nargs=2,default=[18,20])
	parser.add_argument('-f','--figure-font-size',type=int,default=32)
	parser.add_argument('--disable-tex',action='store_true')
	#parser.add_argument('-m','--colormap',default='Oranges')
	parser.add_argument('-m','--colormap',default='BuPu')
	parser.add_argument('--dotcolor',default='yellow')
	parser.add_argument('--dotsize',type=float,default=50)

	return parser.parse_args()

def main():

	# read command line arguments
	args = read_args_from_cmdline()

	expression_file = args.expression_file
	result_file = args.result_file
	go_pickle_file = args.go_pickle_file

	result = pickle.load(open(result_file))
	signatures = result.signatures
	mHG_X_frac = result.mHG_X_frac
	mHG_X_min = result.mHG_X_min
	mHG_L = result.mHG_L

	#go_mHG_X = args.go_mHG_X
	#go_mHG_L = args.go_mHG_L
	#go_pvalue_threshold = args.go_pvalue_threshold
	#go_enrichment_fdr = args.go_enrichment_fdr

	pval_show_below = args.pvalue_show_below
	pval_best = args.pvalue_best
	pval_worst = args.pvalue_worst
	pval_ticks = args.pvalue_ticks
	reverse_order = args.reverse_order

	dim1, dim2 = args.figure_dimensions
	font_size = args.figure_font_size
	cmap = args.colormap
	dotcolor = args.dotcolor
	dotsize = args.dotsize

	max_name_length = args.GO_term_max_name_length
	omit_acc = True
	if args.GO_term_show_id: omit_acc = False

	# read expression data
	print "Reading expression data...", ; sys.stdout.flush()
	genes,samples,E = common.read_expression(expression_file)
	print 'matrix dimensions:', E.shape; sys.stdout.flush()

	# read GO-PCA signatures
	print "Reading GO-PCA signatures...", ; sys.stdout.flush()
	signatures = pickle.load(open(signature_file))
	print "read %d signatures." %(len(signatures)); sys.stdout.flush()

	# read GO pickle
	print "Reading GO data...", ; sys.stdout.flush()
	GO = pickle.load(open(go_pickle_file))
	print "read %d terms, %d annotations." %(len(GO.terms),len(GO.annotations)); sys.stdout.flush()

	# constructing signatures
	print "Constructing signatures...", ;sys.stdout.flush()
	q = len(signatures)
	n = E.shape[1]
	S = np.zeros((q,n),dtype=np.float64)
	labels = []
	for i,sig in enumerate(signatures):
		count = ' (%d)' %(sig.K)
		labels.append(GO.terms[sig.term[0]].get_pretty_format(omit_acc=omit_acc,max_name_length=max_name_length) + count)
		sig = common.get_signature_expression(genes,E,sig.genes)
		S[i,:] = sig
	print 'done!'; sys.stdout.flush()

	# perform hierarchical clustering on signatures (rows)
	print "Clustering of signatures...", ;sys.stdout.flush()
	distxy = squareform(pdist(S, metric='correlation'))
	print 'done!'; sys.stdout.flush()

	# reorder signatures according to clustering
	R = dendrogram(linkage(distxy, method='average'),no_plot=True)
	order = np.int64([int(l) for l in R['ivl']])
	if reverse_order:
		order = order[::-1]
	S = S[order,:]
	labels = [labels[i] for i in order]
	signatures = [signatures[i] for i in order]

	# perform PCA
	print "Performing PCA...", ;sys.stdout.flush()
	test_pc = max(abs(sig.pc) for sig in signatures)
	M_pca = PCA(n_components=test_pc)
	M_pca.fit(E.T)
	print 'done!'; sys.stdout.flush()
	frac = M_pca.explained_variance_ratio_
	cum_frac = np.cumsum(frac)
	print "Cumulative fraction variance explained per PC:"
	print cum_frac
	sys.stdout.flush()

	# test association
	W = M_pca.components_.T
	all_genes = set(genes)
	q = len(signatures)
	A = np.zeros((q,2*test_pc),dtype=np.float64)
	p = E.shape[0]
	matrix = np.zeros((p+1,p+1),dtype=np.longdouble)
	all_genes = set(genes)
	for pc in range(test_pc):
		a_pos = np.argsort(-W[:,pc])
		a_neg = np.argsort(W[:,pc])
		for i,sig in enumerate(signatures):
			term_genes = GO.get_goterm_genes(sig.term[0]) & all_genes
			K = len(term_genes)
		
			v = np.zeros(p,dtype=np.uint8)
			for g in term_genes:
				idx = misc.bisect_index(genes,g)
				v[idx] = 1

			v_sorted = np.ascontiguousarray(v[a_pos])
			threshold,_,pval = mHG_test(v_sorted,p,K,go_mHG_L,go_mHG_X,mat=matrix)
			A[i,pc*2] = -np.log10(pval)
			#if pc == 10 and sig.term[0] == 'GO:0006613':
			#	print -np.log10(pval)
			#	print len(term_genes),np.nonzero(v_sorted)

			v_sorted = np.ascontiguousarray(v[a_neg])
			threshold,_,pval = mHG_test(v_sorted,p,K,go_mHG_L,go_mHG_X,mat=matrix)
			A[i,pc*2+1] = -np.log10(pval)
			#if pc == 10 and sig.term[0] == 'GO:0006613':
			#	print -np.log10(pval)
			#	print len(term_genes),np.nonzero(v_sorted)

	# prepare for plotting
	rc('font',family='serif',size=font_size)
	rc('figure',figsize=(dim1,dim2))
	plt.cla()

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

c="""

#cmap = cm.get_cmap('hot_r')
#cmap = cm.get_cmap('Reds')
#cmap = cm.get_cmap('autumn_r')
#cmap = cm.get_cmap('jet')
#cmap = cm.get_cmap('hot_r')
#cmap.set_bad(color='white')

#H[H<6] = -1.0
#H[H<-9] = -9
#H[H>9] = 9
H[np.absolute(H)<6] = np.nan
plt.imshow(H,interpolation='none',vmin=0,vmax=10,cmap='Oranges',zorder=20)

for p in range(max_pc-1):
    #plt.plot([p*2+0.5,p*2+0.5],[-0.5,m-0.5],color='gray',linewidth=1.0)
    plt.plot([p*2+1.5,p*2+1.5],[-0.5,m-0.5],color='gray',linewidth=1.0,zorder=50)

#plt.gca().spines['right'].set_visible(False)
#plt.gca().spines['bottom'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('top')
    
#plt.plot([-.5,210.5],[0.5,0.5],color='black',linewidth=2)
plt.gca().xaxis.tick_top()
plt.gca().xaxis.set_label_position('top')
plt.xticks(np.arange(0,max_pc*2,2)+0.5,np.arange(max_pc)+1,size='small')
plt.xlabel(r'\textbf{Principal Components}',labelpad=10,size='small')
plt.yticks(np.arange(m),labels,size='x-small')
plt.ylabel(r'\textbf{Signatures}',size='small')
plt.xlim(-0.5,max_pc*2-0.5)
#plt.xlim(-0.5,max_pc*2+1.5)
plt.ylim(m-0.5,-0.5)
plt.grid(which='both',axis='y',zorder=-20) # does not work grid cannot have different z-order than axes

cbar = plt.colorbar(pad=0.02,shrink=0.5)
#ticks = [-8,-6,-4,0,4,6,8]
#ticks = [-8,-6,6,8]
ticks = [6,8,10]
cbar.set_ticks(ticks)
#cbar.ax.tick_params(labelsize='medium')
cbar.set_label(r"$\bm{-\log_{10}}$ p-value")

plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('top')


plt.show()
"""
