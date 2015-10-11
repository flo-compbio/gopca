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
import cPickle as pickle

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

from gopca import common
from genometools import misc

def read_args_from_cmdline():

	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-g','--gopca-file',required=True)
	parser.add_argument('-o','--output-file',required=True)

	parser.add_argument('-fs','--figure-size',type=float,help='in inches',nargs=2,default=[18,18])
	parser.add_argument('-fr','--figure-resolution',type=int,help='in dpi',default=150)
	parser.add_argument('-ff','--figure-font-size',type=int,help='in pt',default=24)
	parser.add_argument('-fm','--figure-font-family',default='serif')
	parser.add_argument('-fc','--figure-colormap',default='RdBu_r')
	parser.add_argument('-fvn','--figure-vmin',type=float,default=-1.0)
	parser.add_argument('-fvx','--figure-vmax',type=float,default=1.0)

	parser.add_argument('-l','--sig-max-name-len',type=int,default=50)
	parser.add_argument('-co','--colorbar-orientation',default='horizontal')
	parser.add_argument('-ca','--colorbar-anchor',type=float,nargs=2,default=(0.96,1.0))
	parser.add_argument('-cs','--colorbar-shrink',type=float,default=0.3)
	parser.add_argument('-cp','--colorbar-pad',type=float,default=0.015)

	parser.add_argument('-t','--use-tex',action='store_true')
	parser.add_argument('-b','--matplotlib-backend',default=None)
	parser.add_argument('-r','--reverse-signature-order',action='store_true')

	return parser.parse_args()

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	result_file = args.gopca_file
	output_file = args.output_file

	# figure size
	fig_size = args.figure_size
	fig_res = args.figure_resolution

	# figure text
	use_tex = args.use_tex
	fig_font_size = args.figure_font_size
	fig_font_family = args.figure_font_family

	# figure heatmap
	fig_vmin = args.figure_vmin
	fig_vmax = args.figure_vmax
	fig_cmap = args.figure_colormap

	# figure colorbar
	cbar_orient = args.colorbar_orientation
	cbar_anchor = args.colorbar_anchor
	cbar_shrink = args.colorbar_shrink
	cbar_pad = args.colorbar_pad

	mpl_backend = args.matplotlib_backend
	reverse_signature_order = args.reverse_signature_order

	sig_max_name_len = args.sig_max_name_len

	# read GO-PCA result
	result = None
	with open(result_file,'rb') as fh:
		result = pickle.load(fh)

	signatures = result.signatures
	# generate labels
	labels = [sig.get_label(include_id=False,max_name_length=sig_max_name_len) for sig in signatures]
	samples = result.samples
	S = result.S

	C = np.corrcoef(S)
	q = S.shape[0]
	assert C.shape[0] == q

	# clustering of rows (signatures)
	order_rows = common.cluster_rows(S,invert=reverse_signature_order)
	C = C[order_rows,:]
	labels = [labels[idx] for idx in order_rows]

	# use the same ordering for columns
	C = C[:,order_rows]

	# plotting
	print 'Plotting...', ; sys.stdout.flush()

	import matplotlib as mpl
	if mpl_backend is not None:
		mpl.use(mpl_backend)
	import matplotlib.pyplot as plt
	from matplotlib import rc

	#if plot_in_notebook:
	#	from IPython import get_ipython
	#	ipython = get_ipython()
	#	ipython.magic('matplotlib inline')

	if use_tex: rc('text',usetex=True)
	rc('font',family=fig_font_family,size=fig_font_size)
	rc('figure',figsize=(fig_size[0],fig_size[1]))
	rc('savefig',dpi=fig_res)

	# plotting
	plt.imshow(C,interpolation='none',aspect='auto',vmin=fig_vmin,vmax=fig_vmax,cmap=fig_cmap)

	plt.xticks(())

	minint = int(fig_vmin)
	maxint = int(fig_vmax)
	cbticks = np.arange(minint,maxint+0.01,0.5)
	cb = plt.colorbar(orientation=cbar_orient,shrink=cbar_shrink,pad=cbar_pad,ticks=cbticks,use_gridspec=False,anchor=cbar_anchor)
	cb.ax.tick_params(labelsize='small')
	cb.set_label('Pearson Correlation',size='small')

	q,n = S.shape
	plt.yticks(np.arange(q),labels,size='x-small')
	plt.ylabel('Signatures')
	plt.xlabel('Signatures')
	print 'done!'; sys.stdout.flush()

	print 'Saving to file...', ; sys.stdout.flush()
	#plt.gcf().set_size_inches(fig_dim[0],fig_dim[1])
	plt.savefig(output_file,bbox_inches='tight')
	print 'done!'; sys.stdout.flush()

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
