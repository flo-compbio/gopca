#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

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

"""

GO-PCA: An Unsupervised Method for Exploring Gene Expression Data \
        Using Prior Knowledge

This script runs GO-PCA and stores the result in a Python "pickle".

Parameters
----------
module_level_variable1 : int
    Module level variables may be documented in either the ``Attributes``
    section of the module docstring, or in an inline docstring immediately
    following the variable.

Example
-------

::

    $ go-pca.py -e [expression_file] -a [annotation_file] -t [ontology-file] -o [output-file]

"""

import sys
import os
import argparse
import logging

import numpy as np

from gopca import common
from gopca.go_pca_objects import GOPCAArgumentParser,GOPCAConfig,GOPCA

def main(args=None):

    # read command line options
    if args is None:
        parser = GOPCAArgumentParser()
        args = parser.parse_args()

    # input files
    expression_file = args.expression_file
    ontology_file = args.ontology_file
    annotation_file = args.annotation_file

    # output file
    output_file = args.output_file

    # log file
    log_file = args.log_file

    # logging parameters
    quiet = args.quiet
    verbose = args.verbose

    log_level = logging.INFO
    if quiet:
        log_level = logging.WARNING
    elif verbose:
        log_level = logging.DEBUG

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

    # for automatically determining the number of PCs
    seed = args.seed
    pc_permutations = args.pc_permutations
    pc_zscore_thresh = args.pc_zscore_thresh

    # intialize logger
    logger = common.get_logger(log_file,log_level)

    ### checks
    assert n_components is None or (isinstance(n_components,int) and n_components >= 0)
    assert isinstance(quiet,bool)
    assert isinstance(verbose,bool)

    # make sure input files exist
    assert os.path.isfile(expression_file)
    assert os.path.isfile(annotation_file)
    if ontology_file is not None:
        assert os.path.isfile(ontology_file)

    # disable global filter if no ontology is provided
    if ontology_file is None:
        logger.warning('Disabling global filter, since no ontology file was provided.')
        disable_global_filter = True

    # generate seed
    if seed is None:
        seed = np.random.randint(int(1e9))
    
    # initialize GO-PCA configuration
    conf_dict = dict([[k,locals()[k]] for k in GOPCAConfig.valid_params])
    config = GOPCAConfig(logger,params=conf_dict)

    # initialize GO-PCA
    M = GOPCA(logger=logger,config=config)

    # read expression data
    M.read_expression(expression_file)

    # filter for most variable genes
    if sel_var_genes > 0:
        M.filter_genes_by_variance(sel_var_genes)

    # estimate the number of PCs (if n_components is set to zero)
    if n_components is None or n_components == 0:
        M.estimate_n_components()
        if M.D == 0:
            logger.error('The estimated number of non-trivial principal components is zero!')
            return 1

    # setting mHG_L to 0 will set the parameter to the default value (= the number of genes / 8)
    if mHG_L == 0:
        L = int(M.p/8.0)
        logger.info('Setting mHG_L to %d.', L)
        M.config.mHG_L = L

    # read ontology
    if ontology_file is not None:
        M.read_ontology(ontology_file)

    # read annotations
    M.read_annotations(annotation_file)

    # run GO-PCA!
    gopca_result = M.run()

    # save output to file
    logger.info('Saving result to file "%s"...', output_file)
    gopca_result.save(output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
