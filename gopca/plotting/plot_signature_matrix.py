import sys
import os
import argparse
import cPickle as pickle

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc

if __name__ == '__main__' and __package__ is None:
	# allow explicit relative imports in executable script
	# source: http://stackoverflow.com/a/6655098
	parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	parent_dir = os.path.dirname(parent_dir)
	print parent_dir; sys.stdout.flush()
	sys.path.insert(1, parent_dir)
	import gopca
	__package__ = 'gopca'

import gopca
from gopca import common
from genometools import misc

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-e','--expression-file',required=True)
	parser.add_argument('-s','--signature-file',required=True)
	parser.add_argument('-g','--go-pickle-file',required=True)

	parser.add_argument('-d','--figure-dimensions',type=float,nargs=2)
	parser.add_argument('-r','--figure-resolution',type=int,help='in dpi',default=150)
	parser.add_argument('-f','--figure-font-size',type=int,default=32)
	parser.add_argument('--vmin',type=float,default=-3.0)
	parser.add_argument('--vmax',type=float,default=3.0)

	parser.add_argument('--disable-sample-clustering',action='store_true')

	return parser.parse_args()

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	expression_file = args.expression_file
	signature_file = args.signature_file
	go_pickle_file = args.go_pickle_file

	fig_dim = args.figure_dimensions
	fig_res = args.figure_resolution
	fig_fontsize = args.figure_font_size

	vmin = args.vmin
	vmax = args.vmax

	# read expression data
	genes,samples,E = common.read_expression(expression_file)

	# read GO pickle
	GO = None
	with open(go_pickle_file) as fh:
		GO = pickle.load(fh)

	# read signatures
	signatures = None
	with open(signature_file) as fh:
		signatures = pickle.load(fh).signatures

	# constructing signatures
	S = np.row_stack([common.get_signature_expression(genes,E,sig.genes) for sig in signatures])
	print S.shape; sys.stdout.flush()
	labels = [common.get_signature_label(GO,sig) for sig in signatures]
	q,n = S.shape
	#print m

	# clustering of rows (signatures)
	print 'Clustering of signatures...', ; sys.stdout.flush()
	distxy = squareform(pdist(S, metric='correlation'))
	R = dendrogram(linkage(distxy, method='average'),no_plot=True)
	order_rows = np.int64([int(l) for l in R['ivl']])
	S = S[order_rows,:]
	labels = [labels[idx] for idx in order_rows]
	print 'done!'; sys.stdout.flush()

	if not args.disable_sample_clustering:
		# clustering of columns (samples)
		print 'Clustering of samples...', ; sys.stdout.flush()
		distxy = squareform(pdist(S.T, metric='euclidean'))
		R = dendrogram(linkage(distxy, method='average'),no_plot=True)
		order_cols = np.int64([int(l) for l in R['ivl']])
		S = S[:,order_cols]
		print 'done!'; sys.stdout.flush()

	# plotting
	print 'Plotting...'; sys.stdout.flush()

	rc('font',family='serif',size=fig_fontsize)
	rc('figure',figsize=(fig_dim[0],fig_dim[1]))
	rc('savefig',dpi=fig_res)

	# plotting
	cmap = 'RdBu_r'
	plt.imshow(S,interpolation='none',aspect='auto',vmin=vmin,vmax=vmax,cmap=cmap)

	minint = int(vmin)
	maxint = int(vmax)
	cbticks = np.arange(minint,maxint+0.01,1.0)
	cbar = plt.colorbar(orientation='horizontal',shrink=0.4,pad=0.02,ticks=cbticks,use_gridspec=False,anchor=(-0.3,1.0))
	cbar.ax.tick_params(labelsize='small')
	cbar.set_label('Standardized Expression',size='small')

	plt.yticks(np.arange(q),labels,size='x-small')
	plt.xlabel('Samples (n=%d)' %(n))
	plt.ylabel('Signatures')
	plt.show()

	plt.cla()
	plt.clf()

	rc('font',family='serif',size=fig_fontsize)
	rc('figure',figsize=(fig_dim[0],fig_dim[1]))
	rc('savefig',dpi=fig_res)

	cmap = 'RdBu_r'
	plt.imshow(np.corrcoef(S),interpolation='none',aspect=1.0,vmin=-1,vmax=1,cmap=cmap)
	cbar = plt.colorbar(shrink=0.4,pad=0.015)
	cbar.ax.tick_params(labelsize='small')
	cbar.set_label('Pearson Correlation')
	plt.yticks(np.arange(q),labels,size='x-small')
	plt.xticks(np.arange(q),(),size='small')
	plt.ylabel('Signatures')
	plt.xlabel('Signatures')

	plt.show()
	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
