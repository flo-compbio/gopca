#!/usr/bin/env python

# Copyright (c) 2015, 2016 Florian Wagner
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

    $ go-pca.py -e [expression_file] -a [annotation_file] -t [ontology-file] \
            -o [output-file]

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
# import os
# import argparse
import textwrap
import logging

import genometools
from genometools.expression import ExpMatrix
from genometools.basic import GeneSetCollection
from genometools.ontology import GeneOntology
from .. import GOPCAParams, GOPCA
from .. import util
from . import arguments


def get_argument_parser():

    prog = 'go-pca.py'
    description = 'Run GO-PCA.'
    parser = arguments.get_argument_parser(prog, description)

    file_mv = arguments.file_mv
    # name_mv = arguments.name_mv
    int_mv = arguments.int_mv
    float_mv = arguments.float_mv
    # str_mv = arguments.str_mv

    # separate config file
    g = parser.add_argument_group('Separate configuration file')
    
    g.add_argument(
        '-c', '--config-file', type=str, metavar=file_mv,
        help=textwrap.dedent("""\
            GO-PCA configuration file. Note: The parameter values specified
            as command line arguments (see below) overwrite the
            corresponding values in the configuration file.""")
    )

    # input and output files
    g = parser.add_argument_group('Input and output files')

    g.add_argument(
        '-e', '--expression-file', type=str, required=True, metavar=file_mv,
        help='Tab-separated text file containing the gene expression matrix.'
    )

    g.add_argument(
        '-s', '--gene-set-file', type=str, required=True, metavar=file_mv,
        help='Tab-separated text file containing the gene sets.'
    )

    g.add_argument(
        '-t', '--gene-ontology-file', type=str, metavar=file_mv, default=None,
        help='OBO file containing the Gene Ontology.'
    )

    g.add_argument(
        '-o', '--output-file', type=str, required=True, metavar=file_mv,
        help='Output pickle file (extension ".pickle" is recommended).'
    )

    # reporting options
    arguments.add_reporting_args(parser)

    # input file hash values
    """
    g = parser.add_argument_group(
        'Input file hash values (can be used for validation purposes)'
    )

    g.add_argument('-he', '--expression-file-hash', metavar = str_mv,
            help = 'MD5 hash for the expression file.')

    g.add_argument('-hs', '--gene-set-file-hash', metavar = str_mv,
            help = 'MD5 hash for the gene set file.')

    g.add_argument('-ht', '--gene-ontology-file-hash', metavar = str_mv,
            help = 'MD5 hash for the gene ontology file.')
    """

    # GO-PCA parameters
    g = parser.add_argument_group('GO-PCA parameters ([] = default value)')

    g.add_argument(
        '-D', '--n-components', type=int, metavar=int_mv, default=-1,
        help=textwrap.dedent("""\
            Number of principal components to test
            (-1 = determine automatically using a permutation test). [-1]"""))

    g.add_argument(
        '-G', '--sel-var-genes', type=int, metavar=int_mv, default=0,
        help=textwrap.dedent("""\
            Variance filter: Keep G most variable genes (0 = off). [0]"""))

    g.add_argument(
        '-P', '--pval-thresh', type=float, metavar=float_mv,
        help=textwrap.dedent("""\
            P-value threshold for GO enrichment test. [%s]
            """ % '%(default).1e'))

    g.add_argument(
        '-E', '--escore-thresh', type=float, metavar=float_mv,
        help=textwrap.dedent("""\
            E-score threshold for GO enrichment test. [%s]
            """ % '%(default).1f'))

    g.add_argument(
        '-R', '--sig-corr-thresh', type=float, metavar=float_mv,
        help=textwrap.dedent("""\
            Correlation threshold used in generating signatures. [%s]
            """ % '%(default).2f'))

    g.add_argument(
        '-Xf', '--mHG-X-frac', type=float, metavar=float_mv,
        help=textwrap.dedent("""\
            X_frac parameter for GO enrichment test. [%s]
            """ % '%(default).2f'))

    g.add_argument(
        '-Xm', '--mHG-X-min', type=int, metavar=int_mv,
        help=textwrap.dedent("""\
            X_min parameter for GO enrichment test. [%s]
            """ % '%(default)d'))

    g.add_argument(
        '-L', '--mHG-L', type=int, metavar=int_mv,
        help=textwrap.dedent("""\
            L parameter for GO enrichment test
            (0 = "off"; -1 = # genes / 8). [%s]""" % '%(default)d'))

    g.add_argument(
        '--escore-pval-thresh', type=float, metavar=float_mv,
        help=textwrap.dedent("""\
            P-value threshold for XL-mHG E-score calculation ("psi"). [%s]
            """ % '%(default).1e'))

    g = parser.add_argument_group('Manually disable the GO-PCA filters')

    g.add_argument(
        '--no-local-filter', action='store_true',
        help='Disable the "local" filter.')

    g.add_argument(
        '--no-global-filter', action='store_true',
        help='Disable the "global" filter (if -t is specified).')

    # legacy options
    g = parser.add_argument_group('Legacy options')

    g.add_argument(
        '--go-part-of-cc-only', action='store_true',
        help='Only propagate "part of" GO relations for the CC domain.')

    # for automatically determining the number of PCs
    # (only used if -D is not specified)
    g = parser.add_argument_group(
        'Parameters for automatically determining the number of PCs '
        'to test ([] = default value)')

    g.add_argument(
        '-ps', '--pc-seed', type=int, metavar=int_mv,
        help=textwrap.dedent("""\
            Random number generator seed (-1 = arbitrary value). [%s]
            """ % '%(default)d'))

    g.add_argument(
        '-pp', '--pc-num-permutations', type=int, metavar=int_mv,
        help='Number of permutations. [%s]' % '%(default)d')

    g.add_argument(
        '-pz', '--pc-zscore-thresh', type=float, metavar=float_mv,
        help='Z-score threshold. [%s]' % '%(default).2f')

    g.add_argument(
        '-pm', '--pc-max-components', type=int, metavar=int_mv,
        help=textwrap.dedent("""\
            Maximum number of PCs to test (0 = no maximum). [%s]
            """ % '%(default)d'))

    # check that the GO-PCA parameter names match the argument names
    # for p in GOPCAParams.param_defaults:
    #    assert p in dir(args)

    # set the argument default values to the parameter defaults stored in
    # the GOPCAParams class
    defaults = GOPCA.get_param_defaults()
    defaults.update(GOPCAParams.get_param_defaults())
    parser.set_defaults(**defaults)

    return parser


def main(args=None):
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
 
    Raises
    ------
    SystemError
        If the version of the Python interpreter is not >= 2.7.
    """
    vinfo = sys.version_info
    if not (vinfo >= (2, 7)):
        raise SystemError('Python interpreter version >= 2.7 required, '
                          'found %d.%d instead.' % (vinfo.major, vinfo.minor))

    if args is None:
        # read arguments from the command line
        parser = get_argument_parser()

        # parse first with default options, in case "--help" is specified
        # ("--help" causes the program to exit at this point)
        args = parser.parse_args()

        # now remove the defaults and parse again
        # (removing the defaults is important so that we know which values
        # were specified by the user)
        no_defaults = dict([p, None] for p in GOPCA.get_param_defaults())
        no_defaults2 = dict([p, None] for p in GOPCAParams.get_param_defaults())
        no_defaults.update(no_defaults2)
        parser.set_defaults(**no_defaults)
        args = parser.parse_args()

    # reporting options
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # test if we can write to log_file?

    # configure root logger
    logger = util.get_logger(log_file=log_file, quiet=quiet)

    # check if required parameters were specified
    passed = True
    if args.expression_file is None:
        logger.error('No expression file specified!')
        passed = False
    if args.gene_set_file is None:
        logger.error('No gene set file specified!')
        passed = False
    if args.output_file is None:
        logger.error('No output file specified!')
        passed = False
    if not passed:
        logger.error('Not all required parameters were specified.')
        return 1

    # generate configuration
    if args.config_file is not None:
        # read parameter values from config file
        params = GOPCAParams.read_ini(args.config_file)
    else:
        # start with default configuration
        params = GOPCAParams()

    # overwrite parameters specified on the command line
    for p in GOPCAParams.get_param_defaults():
        v = getattr(args, p)
        if v is not None:
            logger.debug('Parameter "%s" specified on command line!', p)
            params.set_param(p, v)

    global_params = GOPCA.get_param_defaults()
    for k in list(global_params.keys()):
        v = getattr(args, k)
        if v is not None:
            logger.debug('Parameter "%s" specified on command line!', p)
            global_params[k] = v

    # read expression file
    matrix = ExpMatrix.read_tsv(args.expression_file)
    logger.info('Expression matrix size: ' +
                '(p = %d genes) x (n = %d samples).', matrix.p, matrix.n)

    if args.sel_var_genes > 0:
        # filter genes by variance
        matrix = matrix.filter_variance(args.sel_var_genes)
    
    # read gene set file
    gene_sets = GeneSetCollection.read_tsv(args.gene_set_file)
    print(args.gene_set_file, gene_sets)
    
    # read ontology file (if supplied)
    gene_ontology = None
    if args.gene_ontology_file is not None:
        p_logger = logging.getLogger(genometools.__name__)
        p_logger.setLevel(logging.ERROR)
        gene_ontology = GeneOntology.read_obo(
            args.gene_ontology_file,
            part_of_cc_only=params.go_part_of_cc_only)
        p_logger.setLevel(logging.NOTSET)
        
    M = GOPCA.simple_setup(matrix, params, gene_sets, gene_ontology,
                          verbose=verbose, **global_params)
    run = M.run()

    if run is None:
        logger.error('GO-PCA run failed!')
        return 1

    # write run to pickle file
    logger.info('Storing GO-PCA run in file "%s"...', args.output_file)
    run.write_pickle(args.output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
