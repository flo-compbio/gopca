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

from gopca import common
from gopca.go_pca_objects import GOPCAConfig,GOPCA
#from gopca.printf import printf

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='GO-PCA')

	###
	### Required arguments
	###

	# input files
	parser.add_argument('-e','--expression-file',required=True)
	parser.add_argument('-a','--annotation-file',required=True)
	parser.add_argument('-t','--ontology-file',default=None)

	# output file
	parser.add_argument('-o','--output-file',required=True)

	# number of principal components to test
	parser.add_argument('-D','--principal-components',type=int,required=True)

	###
	### Optional arguments
	###

	# log file
	parser.add_argument('-l','--log-file',default=None)

	# GO-PCA parameters
	parser.add_argument('-P','--pval-thresh',type=float,default=1e-6) # p-value threshold for GO enrichment
	parser.add_argument('-F','--msfe-thresh',type=float,default=2.0) # MSFE (fold enrichment) threshold
	parser.add_argument('-R','--sig-corr-thresh',type=float,default=0.5) # correlation threshold for signature genes

	parser.add_argument('-Xf','--mHG-X-frac',type=float,default=0.25) # 0=off
	parser.add_argument('-Xm','--mHG-X-min',type=int,default=5) # 0=off
	parser.add_argument('-L','--mHG-L',type=int,default=1000) # 0=off
	parser.add_argument('--msfe-pval-thresh',type=float,default=1e-4) # p-value threshold for MSFE calculation

	# variance filter
	parser.add_argument('-N','--select-variable-genes',type=int,default=0)

	# allow filtering to be disabled
	parser.add_argument('--disable-local-filter',action='store_true')
	parser.add_argument('--disable-global-filter',action='store_true')

	### legacy options
	parser.add_argument('--go-part-of-cc-only',action='store_true')

	# output verbosity
	parser.add_argument('-v','--verbosity',type=int,default=3)

	return parser.parse_args()

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	# input files
	expression_file = args.expression_file
	ontology_file = args.ontology_file
	annotation_file = args.annotation_file

	# output file
	output_file = args.output_file

	# log file
	log_file = args.log_file

	# GO-PCA parameters
	n_components = args.principal_components
	pval_thresh = args.pval_thresh
	sig_corr_thresh = args.sig_corr_thresh
	mHG_X_frac = args.mHG_X_frac
	mHG_X_min = args.mHG_X_min
	mHG_L = args.mHG_L
	msfe_pval_thresh = args.msfe_pval_thresh
	msfe_thresh = args.msfe_thresh

	sel_var_genes = args.select_variable_genes

	disable_local_filter = args.disable_local_filter
	disable_global_filter = args.disable_global_filter

	go_part_of_cc_only = args.go_part_of_cc_only
	verbosity = args.verbosity

	# intialize logger
	logger = common.Logger(verbosity,log_file=log_file)

	# make sure input files exist
	assert os.path.isfile(expression_file)
	assert os.path.isfile(annotation_file)
	if ontology_file is not None:
		assert os.path.isfile(ontology_file)

	# checks
	assert isinstance(verbosity,int) and verbosity >= 0

	# disable global filter if no ontology is provided
	if ontology_file is None:
		logger.warning('Disabling global filter, since no ontology file was provided.')
		disable_global_filter = True

	# initialize GO-PCA configuration
	conf_params = ['n_components','pval_thresh','sig_corr_thresh',\
			'mHG_X_frac','mHG_X_min','mHG_L',\
			'msfe_pval_thresh','msfe_thresh',\
			'disable_local_filter','disable_global_filter',\
			'go_part_of_cc_only']
	conf_dict = dict([[k,locals()[k]] for k in conf_params])
	config = GOPCAConfig(logger,**conf_dict)

	# initialize GO-PCA
	M = GOPCA(logger=logger,config=config)

	# read expression data
	M.read_expression(expression_file)

	# setting mHG_L to 0 will "turn off" the effect of the parameter (= set it to the number of genes)
	if mHG_L == 0:
		mHG_L = M.p

	# filter for most variable genes
	if sel_var_genes > 0:
		M.filter_genes_by_variance(sel_var_genes)

	# read ontology
	if ontology_file is not None:
		M.read_ontology(ontology_file)

	# read annotations
	M.read_annotations(annotation_file)

	# run GO-PCA!
	M.run()

	# save output to file
	M.save_result(output_file)

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
