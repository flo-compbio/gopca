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

from gopca import util
from gopca import cli
from gopca import GOPCAConfig, GOPCA

def get_argument_parser():

    prog = 'go-pca.py'
    description = 'Run GO-PCA.'
    parser = cli.get_argument_parser(prog, description)

    file_mv = cli.file_mv
    name_mv = cli.name_mv
    int_mv = cli.int_mv
    float_mv = cli.float_mv
    str_mv = cli.str_mv

    # input and output files
    g = parser.add_argument_group('Input and output files')

    g.add_argument('-e', '--expression-file', required = True,
            metavar = file_mv,
            help = 'Tab-separated text file containing the expression matrix.')

    g.add_argument('-a', '--go-annotation-file', required = True,
            metavar = file_mv,
            help = 'Tab-separated text file containing the GO term ' +
            'annotations.')

    g.add_argument('-t', '--gene-ontology-file', required = False,
            metavar = file_mv,
            help = 'OBO file containing the Gene Ontology.')

    g.add_argument('-o', '--output-file', required = True,
            metavar = file_mv,
            help = 'Output pickle file (extension ".pickle" is recommended).')

    # input file hash values
    g = parser.add_argument_group('Input file hash values')

    g.add_argument('-he', '--expression-file-hash', default = None,
            metavar = str_mv,
            help = 'MD5 hash for the experssion file.')

    g.add_argument('-ha', '--go-annotation-file-hash', default = None,
            metavar = str_mv,
            help = 'MD5 hash for the GO annotation file.')

    g.add_argument('-ht', '--gene-ontology-file-hash', default = None,
            metavar = str_mv,
            help = 'MD5 hash for the gene ontology file.')

    # GO-PCA parameters
    g = parser.add_argument_group('GO-PCA parameters')
    g.add_argument('-D', '--principal-components', type = int, default = 0,
            metavar = int_mv,
            help = 'Number of principal components to test (0 = automatic).')

    g.add_argument('-G', '--select-variable-genes', type = int, default = 0,
            metavar = int_mv,
            help = 'Variance filter: Keep G most variable genes (0 = off).')

    g.add_argument('-P', '--pval-thresh', type = float, default = 1e-6,
            metavar = float_mv,
            help = 'P-value threshold for GO enrichment.')

    g.add_argument('-E', '--escore-thresh', type = float, default = 2.0,
            metavar = float_mv,
            help = 'E-score threshold for GO enrichment.')

    g.add_argument('-R', '--sig-corr-thresh', type = float, default = 0.5,
            metavar = float_mv,
            help = 'Threshold for correlation with seed for signature genes.')

    g.add_argument('-Xf', '--mHG-X-frac', type = float, default = 0.25,
            metavar = float_mv,
            help = 'X_frac parameter for GO enrichment (=> XL-mHG''s X).')

    g.add_argument('-Xm', '--mHG-X-min', type = int, default = 5,
            metavar = int_mv,
            help = 'X_min parameter for GO enrichment (=> XL-mHG''s X).')

    g.add_argument('-L', '--mHG-L', type = int, default = None,
            metavar = int_mv,
            help='L parameter for GO enrichment (0 = off; None = #genes / 8).')

    g.add_argument('--escore-pval-thresh', type = float, default = 1e-4,
            metavar = float_mv,
            help = 'P-value threshold for XL-mHG E-score calculation (= psi).')

    g.add_argument('--no-local-filter', action = 'store_true',
            help = 'Disable GO-PCA''s "local" filter.')

    g.add_argument('--no-global-filter', action = 'store_true',
            help = 'Disable GO-PCA''s "global" filter (if -t is specified).')

    # for automatically determining the number of PCs
    # (only used if -D is not specified)
    g.add_argument('-ps', '--pc-seed', type = int, default = 0,
            metavar = int_mv,
            help = 'Random number generator seed (for determining # PCs).')

    g.add_argument('-pp', '--pc-permutations', type = int, default = 15,
            metavar = int_mv,
            help = 'Number of permutations.')

    g.add_argument('-pz', '--pc-zscore-thresh', type = float, default = 2.0,
            metavar = float_mv,
            help = 'Z-score threshold.')

    # legacy options
    g.add_argument('--go-part-of-cc-only', action = 'store_true',
            help = 'Only propagate "part of" GO relations for the CC domain.')

    # reporting options
    cli.add_reporting_args(parser)

    return parser

def main(args = None):
    """Run GO-PCA and store the result in a `pickle` file.

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
    gene_ontology_file = args.gene_ontology_file
    go_annotation_file = args.go_annotation_file

    # input file hash values
    expression_file_hash = args.expression_file_hash
    gene_ontology_file_hash = args.gene_ontology_file_hash
    go_annotation_file_hash = args.go_annotation_file_hash

    # output file
    output_file = args.output_file

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
    no_local_filter = args.no_local_filter
    no_global_filter = args.no_global_filter
    go_part_of_cc_only = args.go_part_of_cc_only

    # for automatically determining the number of PCs
    pc_seed = args.pc_seed
    pc_permutations = args.pc_permutations
    pc_zscore_thresh = args.pc_zscore_thresh

    # reporting options
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # test if we can write to log_file?

    # configure root logger
    logger = util.get_logger(log_file = log_file, quiet = quiet,
            verbose = verbose)

    # create GO-PCA configuration
    params = {}
    for p in sorted(GOPCAConfig.param_names):
        params[p] = locals()[p]

    # why is the line below not working?
    #params = dict([k, locals()[k]] for k in sorted(GOPCAConfig.param_names))

    config = GOPCAConfig(params)

    M = GOPCA(config)
    G = M.run()
    G.save(output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
