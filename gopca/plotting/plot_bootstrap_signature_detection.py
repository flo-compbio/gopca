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

from genometools import misc
from goparser import GOParser
from gopca import common

def read_args_from_cmdline():

	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-g','--gopca-file',required=True)
	parser.add_argument('-b','--bootstrap-gopca-file',required=True)
	parser.add_argument('-t','--ontology-file',required=True)
	parser.add_argument('-o','--output-file',required=True)

	parser.add_argument('--part-of-cc-only',action='store_true')

	parser.add_argument('-fs','--figure-size',type=float,help='in inches',nargs=2,default=[18,12])
	parser.add_argument('-fr','--figure-resolution',type=int,help='in dpi',default=150)
	parser.add_argument('-ffs','--figure-font-size',type=int,help='in pt',default=32)
	parser.add_argument('-fff','--figure-font-family',default='serif')
	parser.add_argument('-ft','--figure-use-tex',action='store_true')
	parser.add_argument('-fb','--figure-backend',default=None)
	parser.add_argument('-fy','--figure-ylim',type=float,default=None)

	return parser.parse_args()

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	gopca_file = args.gopca_file
	bootstrap_gopca_file = args.bootstrap_gopca_file
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

	fig_ylim = args.figure_ylim

	# other
	mpl_backend = args.figure_backend

	gopca_result = common.read_gopca_result(gopca_file)
	samples = gopca_result.samples
	signatures = gopca_result.signatures
	S = gopca_result.S
	q = len(signatures)

	GO = GOParser()
	GO.parse_ontology(ontology_file,part_of_cc_only)

	# determine set of "related" GO terms for each signature
	sig_term_ids = [sig.term[0] for sig in signatures]
	related_term_ids = [set([id_]) | GO.terms[id_].ancestors | GO.terms[id_].children for id_ in sig_term_ids]

	bootstrap_result = common.read_gopca_result(bootstrap_gopca_file)
	T = bootstrap_result.T
	R_ident = np.zeros((T,q),dtype=np.float64)
	R_rel = np.zeros((T,q),dtype=np.float64)

	for j,result in enumerate(bootstrap_result.gopca_results):
		for i,sig in enumerate(result.signatures):
			for i_ref in range(q):
				if sig.term[0] in related_term_ids[i_ref]:
					R_rel[j,i_ref] = 1.0
					if sig.term[0] == sig_term_ids[i_ref]:
						R_ident[j,i_ref] = 1.0
	R_ident = 100*np.mean(R_ident,axis=0)
	R_rel = 100*np.mean(R_rel,axis=0)

	bin_edges = np.r_[np.arange(0,100.0,25.0),100.1]
	ident_binned = np.bincount(np.digitize(R_ident,bin_edges)-1)
	rel_binned = np.bincount(np.digitize(R_rel,bin_edges)-1)

	import matplotlib as mpl
	if mpl_backend is not None:
		mpl.use(mpl_backend)
	import matplotlib.pyplot as plt
	from matplotlib import rc

	if fig_use_tex: rc('text',usetex=True)
	rc('font',family=fig_font_family,size=fig_font_size)
	rc('figure',figsize=(fig_size[0],fig_size[1]))
	rc('savefig',dpi=fig_res)

	b1 = plt.bar(np.arange(bin_edges.size-1)-0.4,ident_binned,width=0.38,edgecolor='none',color='skyblue')
	b2 = plt.bar(np.arange(bin_edges.size-1)+0.02,rel_binned,width=0.38,edgecolor='none',color='gold')
	xticklabels = [r'%.0f - %.0f' %(bin_edges[i],bin_edges[i+1]) for i in range(bin_edges.size-1)]
	plt.xticks(np.arange(bin_edges.size-1),xticklabels)
	plt.tick_params(axis='x',length=0)
	if fig_use_tex:
		plt.xlabel(r'Detection Rate (\% Bootstrap Samples)')
		plt.ylabel(r'\# Signatures')
	else:
		plt.xlabel(r'Detection Rate (% Bootstrap Samples)')
		plt.ylabel('# Signatures')
	plt.legend([b1,b2],['Identical GO Term', 'Related GO Term'],loc='upper right',fontsize='small')
	if fig_ylim is not None:
		plt.ylim(0,fig_ylim)

	common.simpleaxis(plt.gca())
	print 'Saving to file...', ; sys.stdout.flush()
	plt.savefig(output_file,bbox_inches='tight')
	print 'done!'; sys.stdout.flush()

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
