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

"""This script runs GO-PCA and stores the result in a Python "pickle".

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

from genometools import misc
from gopca import GOPCAInput,GOPCA
#from gopca.go_pca_objects import GOPCAArgumentParser,GOPCAConfig,GOPCA

def get_argument_parser():

    description = 'GO-PCA: An Unsupervised Method to Explore Gene Expression Data Using Prior Knowledge'
    formatter_class = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(description = description,
            formatter_class = formatter_class, add_help = False)


    helpgroup = parser.add_argument_group('Help')
    helpgroup.add_argument('-h','--help',action='help',help='show this help message and exit')

    iogroup = parser.add_argument_group('Input and output files')
    # input files
    iogroup.add_argument('-e','--expression-file',required=True,
            help='Tab-separated text file containing the expression matrix.')
    iogroup.add_argument('-a','--go-annotation-file',required=True,
            help='Tab-separated text file containing the GO term annotations.')
    iogroup.add_argument('-t','--ontology-file',required=False,
            help='OBO file containing the Gene Ontology.')
    # output file
    iogroup.add_argument('-o','--output-file',required=True,
            help='Output pickle file (extension ".pickle" is recommended).')

    # reporting
    repgroup = parser.add_argument_group('Reporting options')
    repgroup.add_argument('-l','--log-file',default=None,
            help='Path of log file (if specified, report to stdout AND file.')
    repgroup.add_argument('-q','--quiet',action='store_true',
            help='Only output errors and warnings.')
    repgroup.add_argument('-v','--verbose',action='store_true',
            help='Enable verbose output. Ignored if --quiet is specified.')

    # GO-PCA parameters
    gg = parser.add_argument_group('GO-PCA parameters')
    gg.add_argument('-D','--principal-components',type=int, default=0,
            help='Number of principal components to test (0 = automatic).')

    gg.add_argument('-G','--select-variable-genes',type=int,default=0,\
            help='Variance filter: Keep G most variable genes (0 = off).')

    gg.add_argument('-P','--pval-thresh',type=float,default=1e-6,
            help='P-value threshold for GO enrichment.')

    gg.add_argument('-E','--escore-thresh',type=float,default=2.0,
            help='E-score threshold for GO enrichment.')

    gg.add_argument('-R','--sig-corr-thresh',type=float,default=0.5,
            help='Threshold for correlation with seed for signature genes.')

    gg.add_argument('-Xf','--mHG-X-frac',type=float,default=0.25,
            help='X_frac parameter for GO enrichment (=> XL-mHG''s X).')

    gg.add_argument('-Xm','--mHG-X-min',type=int,default=5,
            help='X_min parameter for GO enrichment (=> XL-mHG''s X).')

    gg.add_argument('-L','--mHG-L',type=int,default=None,
            help='L parameter for GO enrichment (0 = off; None = #genes/8).')

    gg.add_argument('--escore-pval-thresh',type=float,default=1e-4,
            help='P-value threshold for XL-mHG E-score calculation (= psi).')

    # GO term filters
    gg.add_argument('--disable-local-filter',action='store_true',
            help='Disable GO-PCA''s "local" filter.')
    gg.add_argument('--disable-global-filter',action='store_true',
            help='Disable GO-PCA''s "global" filter (if -t is specified).')

    # for automatically determining the number of PCs
    # (only used if -D is unset)
    gg.add_argument('-s','--seed',type=int,default=None,
            help='Random number generator seed (None = not fixed).')
    gg.add_argument('-pp','--pc-permutations',type=int,default=15,
            help='Number of permutations.')
    gg.add_argument('-pz','--pc-zscore-thresh',type=float,default=2.0,
            help='Z-score threshold.')

    # legacy options
    gg.add_argument('--go-part-of-cc-only',action='store_true',
            help='Only propagate "part of" GO relations for the CC domain.')

    return parser

def main(args=None):
    """Run GO-PCA and store the result in a Python "pickle".

    Parameters
    ----------
    args: argparse.Namespace object, optional
        The argument values. If not specified, the values will be obtained by
        parsing the command line arguments using the `argparse` module.

    Returns
    -------
    int
        Exit code (0 if no error occurred).
 
    """
    # read command line options
    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    # input files
    expression_file = args.expression_file
    ontology_file = args.ontology_file
    go_annotation_file = args.go_annotation_file

    # output file
    output_file = args.output_file

    # logging parameters
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

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

    # configure root logger
    log_level = logging.INFO
    if quiet:
        log_level = logging.WARNING
    elif verbose:
        log_level = logging.DEBUG

    logger = misc.configure_logger('', log_file = log_file,
            log_level = log_level)

    # generate random seed (if not provided)
    if seed is None:
        seed = np.random.randint(int(1e9))

    # collect and validate GO-PCA input
    inpt = GOPCAInput()
    params = dict([[k,locals()[k]] for k in sorted(inpt.param_names)])
    inpt.set_params(params)

    M = GOPCA(inpt)

    output = M.run()

    output.save(output_file)

    """
    # save output to file
    logger.info('Saving result to file "%s"...', output_file)
    gopca_result.save(output_file)
    """

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
