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
from scipy.stats import pearsonr

from genometools import misc
from goparser import GOParser
from gopca import common

def read_args_from_cmdline():

	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-g','--gopca-file',required=True)
	parser.add_argument('-b','--bootstrap-gopca-files',required=True,nargs='+')
	parser.add_argument('-t','--ontology-file',required=True)
	parser.add_argument('-o','--output-file',required=True)

	parser.add_argument('--part-of-cc-only',action='store_true')

	parser.add_argument('-fs','--figure-size',type=float,help='in inches',nargs=2,default=[18,12])
	parser.add_argument('-fr','--figure-resolution',type=int,help='in dpi',default=150)
	parser.add_argument('-ff','--figure-font-size',type=int,help='in pt',default=32)
	parser.add_argument('-fm','--figure-font-family',default='serif')
	parser.add_argument('-ft','--figure-use-tex',action='store_true')
	parser.add_argument('-fb','--figure-backend',default=None)
	#parser.add_argument('-fc','--figure-colormap',default='RdBu_r')
	#parser.add_argument('-fvn','--figure-vmin',type=float,default=-3.0)
	#parser.add_argument('-fvx','--figure-vmax',type=float,default=3.0)

	return parser.parse_args()

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	gopca_file = args.gopca_file
	bootstrap_gopca_files = args.bootstrap_gopca_files
	ontology_file = args.ontology_file
	output_file = args.output_file
	part_of_cc_only = args.part_of_cc_only

	# figure size
	fig_size = args.figure_size
	fig_res = args.figure_resolution

	# figure text
	fig_use_tex = args.figure_use_tex
	fig_font_size = args.figure_font_size
	fig_font_family = args.figure_font_family

	mpl_backend = args.figure_backend

	gopca_result = common.read_gopca_result(gopca_file)
	samples = gopca_result.samples
	signatures = gopca_result.signatures
	S = gopca_result.S

	# order doesn't matter here
	#a = common.cluster_signatures(S)
	# if invert:
	# 	a = a[::-1]
	#signatures = [signatures[i] for i in a]
	#S = S[a,:]

	GO = GOParser()
	GO.parse_ontology(ontology_file,part_of_cc_only)

	# determine set of "acceptable" GO terms for each signature
	sig_term_ids = [sig.term[0] for sig in signatures]
	valid_term_ids = [set([id_]) | GO.terms[id_].ancestors | GO.terms[id_].children for id_ in sig_term_ids]
	#valid_term_ids = [set([id_]) for id_ in sig_term_ids]

	import matplotlib as mpl
	if mpl_backend is not None:
		mpl.use(mpl_backend)
	import matplotlib.pyplot as plt
	from matplotlib import rc

	if fig_use_tex: rc('text',usetex=True)
	rc('font',family=fig_font_family,size=fig_font_size)
	rc('figure',figsize=(fig_size[0],fig_size[1]))
	rc('savefig',dpi=fig_res)

	q = len(signatures)
	k = len(bootstrap_gopca_files)

	repeats = np.zeros(k,dtype=np.float64)
	resample_sizes = np.zeros(k,dtype=np.float64)

	det_med = np.zeros(k,dtype=np.float64)
	det_lq = np.zeros(k,dtype=np.float64)
	det_uq = np.zeros(k,dtype=np.float64)

	corr_med = np.zeros(k,dtype=np.float64)
	corr_lq = np.zeros(k,dtype=np.float64)
	corr_uq = np.zeros(k,dtype=np.float64)

	for h,fn in enumerate(bootstrap_gopca_files):
		# for each bootstrap run
		# we're only interested in averages across bootstrap samples here

		bootstrap_result = common.read_gopca_result(fn)
		resample_sizes[h] = bootstrap_result.resample_size
		T = bootstrap_result.T
		found = np.zeros((T,q),dtype=np.float64)
		maxcorr = np.zeros((T,q),dtype=np.float64) - 1.0
		for j,result in enumerate(bootstrap_result.gopca_results):
			# for each bootstrap sample
			indices = np.int64([samples.index(s) for s in result.samples])
			S_sub = S[:,indices]
			found_repeat = np.zeros(q,dtype=np.float64)
			for i,sig in enumerate(result.signatures):
				c = -1.0
				for i_ref in range(q):
					if sig.term[0] in valid_term_ids[i_ref]:
						found_repeat[i_ref] = 1.0
					r,_ = pearsonr(S_sub[i_ref,:],result.S[i,:])
					maxcorr[j,i_ref] = max(maxcorr[j,i_ref],r)
			found[j,:] = found_repeat
		found = 100*np.mean(found,axis=0)
		maxcorr = np.mean(maxcorr,axis=0)

		det_med[h] = np.median(found)
		det_lq[h] = np.percentile(found,25.0)
		det_uq[h] = np.percentile(found,75.0)
		corr_med[h] = np.median(maxcorr)
		corr_lq[h] = np.percentile(maxcorr,25.0)
		corr_uq[h] = np.percentile(maxcorr,75.0)

	error_kw = {'ecolor': 'orange', 'lw': 3.0, 'capsize': 10, 'mew': 3.0, 'zorder': 50}
	plt.errorbar(np.arange(k),det_med,yerr=[det_med-det_lq,det_uq-det_med],\
			color='orange',**error_kw)
	plt.ylim(0,100.0)
	if fig_use_tex:
		plt.ylabel(r'\parbox{\textwidth}{\centering GO Term Recovery\\(\% Bootstrap Samples)}')
	else:
		plt.ylabel('GO Term Recovery\n(% Bootstrap Samples)')
	common.simpleaxis(plt.gca())

	if fig_use_tex:
		plt.xlabel(r'Sample Size (\% Original Sample Size)')
	else:
		plt.xlabel('Sample Size (% Original Sample Size)')

	ax2 = plt.twinx()
	plt.sca(ax2)

	error_kw = {'ecolor': 'purple', 'lw': 3.0, 'capsize': 10, 'mew': 3.0, 'zorder': 50}
	plt.errorbar(np.arange(k),corr_med,yerr=[corr_med-corr_lq,corr_uq-corr_med],\
			color='purple',**error_kw)

	plt.xticks(np.arange(k),['%.0f' %(s) for s in resample_sizes])
	plt.xlim(-0.5,k-0.5)
	if fig_use_tex:
		plt.ylabel(r'\parbox{\textwidth}{\centering Signature Recovery\\(Max. Correlation)}',color='purple')
	else:
		plt.ylabel('Signature Recovery\n(Max. Correlation)',color='purple')
	plt.yticks(color='purple')
	plt.ylim(0,1.0)
	plt.gca().spines['top'].set_visible(False)

	#for r in range(T):
	#	comps[r] = np.median(C[r,:])
	#	lq[h] = comps[h] - np.percentile(C[h,:],25.0)
	#	uq[h] = np.percentile(C[h,:],75.0) - comps[h]
	print 'Saving to file...', ; sys.stdout.flush()
	#plt.gcf().set_size_inches(fig_dim[0],fig_dim[1])
	plt.savefig(output_file,bbox_inches='tight')
	print 'done!'; sys.stdout.flush()
	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
