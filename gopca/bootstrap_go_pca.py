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
from gopca.bootstrap_go_pca_objects import BootstrapGOPCAArgumentParser,BootstrapGOPCAConfig,BootstrapGOPCA

def main(args=None):

	if args is None:
		parser = BootstrapGOPCAArgumentParser()
		args = parser.parse_args()

	# input files
	expression_file = args.expression_file
	ontology_file = args.ontology_file
	annotation_file = args.annotation_file

	# output file
	output_file = args.output_file

	# log file
	log_file = args.log_file

	# GO-PCA parameters
	sel_var_genes = args.select_variable_genes
	n_components = args.principal_components
	pval_thresh = args.pval_thresh
	sig_corr_thresh = args.sig_corr_thresh
	mHG_X_frac = args.mHG_X_frac
	mHG_X_min = args.mHG_X_min
	mHG_L = args.mHG_L
	escore_pval_thresh = args.escore_pval_thresh
	escore_thresh = args.escore_thresh
	disable_local_filter = args.disable_local_filter
	disable_global_filter = args.disable_global_filter
	go_part_of_cc_only = args.go_part_of_cc_only

	seed = args.seed

	# bootstrap parameters
	repeats = args.repeats
	resample_size = args.resample_size

	# parameters for estimating the number of PCs
	pc_zscore_thresh = args.pc_zscore_thresh
	pc_permutations = args.pc_permutations

	# verbosity
	verbosity = args.verbosity

	# intialize logger
	logger = common.Logger(verbosity,log_file=log_file)

	if seed is None:
		seed = np.random.randint(int(1e9))
	logger.message('Using seed: %d' %(seed))

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
	conf_dict = dict([[k,locals()[k]] for k in BootstrapGOPCAConfig.all_valid_params])
	config = BootstrapGOPCAConfig(logger,**conf_dict)

	# initialize Bootstrap-GO-PCA
	M = BootstrapGOPCA(logger,config)

	# read expression data
	M.read_expression(expression_file)

	# filter for most variable genes
	if sel_var_genes > 0:
		M.filter_genes_by_variance(sel_var_genes)

	# setting mHG_L to 0 will set the parameter to the default value (= the number of genes / 8)
	if mHG_L == 0:
		L = int(M.p/8.0)
		logger.message('Setting mHG_L to %d.' %(L))
		M.config.mHG_L = L


	# read ontology
	if ontology_file is not None:
		M.read_ontology(ontology_file)

	# read annotations
	M.read_annotations(annotation_file)

	# run GO-PCA!
	result = M.run()

	# save output to file
	result.save(output_file)

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
