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
	parser.add_argument('-c','--components',type=int,nargs=2,default=[1,2])
	parser.add_argument('-s','--signature-name',required=True)
	parser.add_argument('-o','--output-file',required=True)


	parser.add_argument('-fs','--figure-size',type=float,help='in inches',nargs=2,default=[16,14])
	parser.add_argument('-fr','--figure-resolution',type=int,help='in dpi',default=150)
	parser.add_argument('-ffs','--figure-font-size',type=int,help='in pt',default=32)
	parser.add_argument('-ffm','--figure-font-family',default='serif')
	parser.add_argument('-fc','--figure-colormap',default='RdBu_r')
	parser.add_argument('-fvn','--figure-vmin',type=float,default=-3.0)
	parser.add_argument('-fvx','--figure-vmax',type=float,default=3.0)
	parser.add_argument('-ft','--figure-use-tex',action='store_true')
	parser.add_argument('-fms','--figure-marker-size',type=float,default=200.0)
	parser.add_argument('-fma','--figure-marker-alpha',type=float,default=0.7)
	#parser.add_argument('-fp','--figure-title-pos',type=float,default=0.95)

	parser.add_argument('-co','--figure-colorbar-orientation',default='horizontal')
	parser.add_argument('-ca','--figure-colorbar-anchor',type=float,nargs=2,default=(0.96,1.0))
	parser.add_argument('-cs','--figure-colorbar-shrink',type=float,default=0.3)
	parser.add_argument('-cp','--figure-colorbar-pad',type=float,default=0.015)

	parser.add_argument('-b','--matplotlib-backend',default=None)

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
	comps = args.components
	sig_name = args.signature_name
	output_file = args.output_file

	# figure size
	fig_size = args.figure_size
	fig_res = args.figure_resolution

	# marker
	fig_marker_size = args.figure_marker_size
	fig_marker_alpha = args.figure_marker_alpha

	# figure text
	fig_use_tex = args.figure_use_tex
	fig_font_size = args.figure_font_size
	fig_font_family = args.figure_font_family

	# figure heatmap
	fig_vmin = args.figure_vmin
	fig_vmax = args.figure_vmax
	fig_cmap = args.figure_colormap

	# figure colorbar
	fig_cbar_orient = args.figure_colorbar_orientation
	fig_cbar_anchor = args.figure_colorbar_anchor
	fig_cbar_shrink = args.figure_colorbar_shrink
	fig_cbar_pad = args.figure_colorbar_pad

	# specific
	#fig_title_pos = args.figure_title_pos
	#fig_subgrid_ratio = args.figure_subgrid_ratio

	mpl_backend = args.matplotlib_backend

	# read GO-PCA result
	result = None
	with open(result_file,'rb') as fh:
		result = pickle.load(fh)

	# find signature selected
	signatures = result.signatures
	term_ids = set([sig.term[0] for sig in signatures])
	sig = None
	if sig_name in term_ids:
		sig = [s for s in signatures if s.term[0] == sig_name]
		assert len(sig) == 1
		sig = sig[0]
	else:
		sig_name = sig_name.lower()
		sig = [s for s in signatures if s.term[3].lower().startswith(sig_name)]
		if len(sig) == 0:
			print >> sys.stderr, 'Error: signature name not found.'
			return 1
		elif len(sig) > 1:
			print >> sys.stderr, 'Error: signature name not unique, matched: %s' %(', '.join([s.term[3] for s in sig]))
			return 1
		sig = sig[0]

	Y = result.Y
	y1 = Y[:,comps[0]-1]
	y2 = Y[:,comps[1]-1]

	# get signature gene expression matrix and cluster rows
	expr = np.mean(common.get_standardized_matrix(sig.E),axis=0)

	# plotting
	import matplotlib as mpl
	if mpl_backend is not None:
		mpl.use(mpl_backend)
	import matplotlib.pyplot as plt
	from matplotlib import rc

	if fig_use_tex: rc('text',usetex=True)
	rc('font',family=fig_font_family,size=fig_font_size)
	rc('figure',figsize=(fig_size[0],fig_size[1]))
	rc('savefig',dpi=fig_res)

	print expr.shape,expr.dtype
	print y1.shape,y2.shape
	plt.scatter(y1,y2,c=expr,s=fig_marker_size,alpha=fig_marker_alpha,cmap=fig_cmap,vmin=fig_vmin,vmax=fig_vmax)

	plt.xticks(())
	plt.xlabel('PC %d Score' %(comps[0]))
	plt.yticks(())
	plt.ylabel('PC %d Score' %(comps[1]))
	plt.title(sig.get_label(include_id=False))

	minint = int(fig_vmin)
	maxint = int(fig_vmax)
	cbticks = np.arange(int(minint),int(maxint)+0.01,1.0)
	cb = plt.colorbar(orientation=fig_cbar_orient,shrink=fig_cbar_shrink,pad=fig_cbar_pad,ticks=cbticks,use_gridspec=False,anchor=fig_cbar_anchor)
	cb.ax.tick_params(labelsize='small')
	cb.set_label('Standardized Expression',size='small')
	simpleaxis(plt.gca())
	#plt.suptitle(sig_label,va='top',y=fig_title_pos)

	#plt.tight_layout()

	print 'Saving to file...', ; sys.stdout.flush()
	plt.savefig(output_file,bbox_inches='tight')
	print 'done!'; sys.stdout.flush()

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
