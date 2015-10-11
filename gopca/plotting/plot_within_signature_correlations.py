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

from gopca import common
from genometools import misc

def read_args_from_cmdline():

	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-g','--gopca-file',required=True)
	parser.add_argument('-o','--output-file',required=True)

	parser.add_argument('-s','--figure-size',type=float,help='in inches',nargs=2,default=[18,18])
	parser.add_argument('-r','--figure-resolution',type=int,help='in dpi',default=150)
	parser.add_argument('-f','--figure-font-size',type=int,help='in pt',default=32)
	parser.add_argument('-m','--figure-font-family',default='serif')
	#parser.add_argument('-xn','--figure-xmin',type=float,default=0.8)
	#parser.add_argument('-xx','--figure-xmax',type=float,default=1.0)

	parser.add_argument('-l','--sig-max-name-len',type=int,default=50)

	parser.add_argument('-t','--use-tex',action='store_true')
	parser.add_argument('-b','--matplotlib-backend',default=None)
	parser.add_argument('-i','--invert-signature-order',action='store_true')

	return parser.parse_args()

def simpleaxis(ax):
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()

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

	# scale
	#fig_xmin = args.figure_xmin
	#fig_xmax = args.figure_xmax

	mpl_backend = args.matplotlib_backend
	invert_signature_order = args.invert_signature_order

	sig_max_name_len = args.sig_max_name_len

	# read GO-PCA result
	result = None
	with open(result_file,'rb') as fh:
		result = pickle.load(fh)

	# generate labels
	signatures = result.signatures
	labels = [sig.get_label(include_id=False,max_name_length=sig_max_name_len) for sig in signatures]
	samples = result.samples
	S = result.S

	# clustering of rows (signatures)
	order_rows = common.cluster_signatures(S,invert=invert_signature_order)
	S = S[order_rows,:]
	labels = [labels[idx] for idx in order_rows]

	groups = []
	for j,sig in enumerate(signatures):
		C = np.corrcoef(sig.E)
		sel = np.triu_indices(C.shape[0],k=1)
		groups.append(C[sel])

	print len(groups),len(labels)

	# plotting
	import matplotlib as mpl
	if mpl_backend is not None:
		mpl.use(mpl_backend)
	import matplotlib.pyplot as plt
	from matplotlib import rc

	if use_tex: rc('text',usetex=True)
	rc('font',family=fig_font_family,size=fig_font_size)
	rc('figure',figsize=(fig_size[0],fig_size[1]))
	rc('savefig',dpi=fig_res)

	boxprops = {'ec':'none', 'fc':'skyblue', 'lw':0}
	medianprops = {'color':'black', 'lw':2.0}
	whiskerprops = {'color':'black'}
	#flierprops = {'color':'gray'}
	flierprops = {'color':'skyblue','markeredgewidth':2.0}

	bp = plt.boxplot(groups,vert=False,labels=labels,flierprops=flierprops,medianprops=medianprops,whiskerprops=whiskerprops,\
           patch_artist=True,boxprops=boxprops)

	plt.setp(bp['boxes'], zorder=100)
	plt.setp(bp['whiskers'], zorder=100)
	plt.setp(bp['fliers'], zorder=100)
	plt.setp(bp['medians'], zorder=200)
	plt.setp(bp['caps'], zorder=100)


	plt.yticks(size='x-small')
	plt.xlabel('Pairwise Correlation of Signature Genes')
	plt.ylabel('Signatures')

	#plt.grid(b=True,which='major',axis='x',lw=3.0,color='red')
	q = S.shape[0]
	plt.grid(b=True,which='major',axis='y',lw=3.0,color='darkgray')
	for x in [-0.5,0,0.5]:
		plt.plot([x,x],[0.5,q+0.5],'-',color='pink',lw=2.0)

	plt.xlim(-1,1)
	#plt.ylim(0.5,q+0.5)
	plt.ylim(q+0.5,0.5)
	simpleaxis(plt.gca())

	print 'Saving to file...', ; sys.stdout.flush()
	#plt.gcf().set_size_inches(fig_dim[0],fig_dim[1])
	plt.savefig(output_file,bbox_inches='tight',zorder=0)
	print 'done!'; sys.stdout.flush()

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
