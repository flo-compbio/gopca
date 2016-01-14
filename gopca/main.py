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
import textwrap

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
    str_type = cli.str_type

    # separate config file
    g = parser.add_argument_group('Separate configuration file')
    
    g.add_argument('-c', '--config-file', type = str_type, metavar = file_mv,
            help = textwrap.dedent("""\
                GO-PCA configuration file. Note: The parameter values specified
                as command line arguments (see below) overwrite the
                corresponding values in the configuration file."""))

    # input and output files
    g = parser.add_argument_group('Input and output files')

    g.add_argument('-e', '--expression-file', type = str_type,
            metavar = file_mv, help = textwrap.dedent("""\
                Tab-separated text file containing the gene expression matrix.
                """))

    g.add_argument('-s', '--gene-set-file', type = str_type,
            metavar = file_mv, help = textwrap.dedent("""\
                Tab-separated text file containing the gene sets.
                """))

    g.add_argument('-t', '--gene-ontology-file', required = False,
            metavar = file_mv, type = str_type,
            help = 'OBO file containing the Gene Ontology.')

    g.add_argument('-o', '--output-file', type = str_type, metavar = file_mv,
            help = 'Output pickle file (extension ".pickle" is recommended).')

    # input file hash values
    g = parser.add_argument_group('Input file hash values (can be used for validation purposes)')

    g.add_argument('-he', '--expression-file-hash', metavar = str_mv,
            help = 'MD5 hash for the experssion file.')

    g.add_argument('-hs', '--gene-set-file-hash', metavar = str_mv,
            help = 'MD5 hash for the gene set file.')

    g.add_argument('-ht', '--gene-ontology-file-hash', metavar = str_mv,
            help = 'MD5 hash for the gene ontology file.')

    # GO-PCA parameters
    g = parser.add_argument_group('GO-PCA parameters ([] = default value)')
    g.add_argument('-D', '--n-components', type = int,
            metavar = int_mv, help = textwrap.dedent("""\
                Number of principal components to test
                (-1 = determine automatically using a permutation test). [%s]
                """ %('%(default)d')))

    g.add_argument('-G', '--sel-var-genes', type = int,
            metavar = int_mv, help = textwrap.dedent("""\
                Variance filter: Keep G most variable genes (0 = off). [%s]
                """ %('%(default)d')))

    g.add_argument('-P', '--pval-thresh', type = float,
            metavar = float_mv, help = textwrap.dedent("""\
                P-value threshold for GO enrichment test. [%s]
                """ %('%(default).1e')))

    g.add_argument('-E', '--escore-thresh', type = float,
            metavar = float_mv, help = textwrap.dedent("""\
                E-score threshold for GO enrichment test. [%s]
                """ %('%(default).1f')))

    g.add_argument('-R', '--sig-corr-thresh', type = float,
            metavar = float_mv, help = textwrap.dedent("""\
                Correlation threshold used in generating signatures. [%s]
                """ %('%(default).2f')))

    g.add_argument('-Xf', '--mHG-X-frac', type = float,
            metavar = float_mv, help = textwrap.dedent("""\
                X_frac parameter for GO enrichment test. [%s]
                """ %('%(default).2f')))

    g.add_argument('-Xm', '--mHG-X-min', type = int,
            metavar = int_mv, help = textwrap.dedent("""\
                X_min parameter for GO enrichment test. [%s]
                """ %('%(default)d')))

    g.add_argument('-L', '--mHG-L', type = int,
            metavar = int_mv, help = textwrap.dedent("""\
                L parameter for GO enrichment test
                (0 = "off"; -1 = # genes / 8). [%s]""" %('%(default)d')))

    g.add_argument('--escore-pval-thresh', type = float,
            metavar = float_mv, help = textwrap.dedent("""\
                P-value threshold for XL-mHG E-score calculation ("psi"). [%s]
                """ %('%(default).1e')))

    g = parser.add_argument_group('Manually disable the GO-PCA filters')
    g.add_argument('--no-local-filter', action = 'store_true',
            help = 'Disable the "local" filter.')

    g.add_argument('--no-global-filter', action = 'store_true',
            help = 'Disable the "global" filter (if -t is specified).')

    # legacy options
    g = parser.add_argument_group('Legacy options')

    g.add_argument('--go-part-of-cc-only', action = 'store_true',
            help = 'Only propagate "part of" GO relations for the CC domain.')


    # for automatically determining the number of PCs
    # (only used if -D is not specified)
    g = parser.add_argument_group('Parameters for automatically determining the number of PCs to test ([] = default value)')
    g.add_argument('-ps', '--pc-seed', type = int,
            metavar = int_mv, help = textwrap.dedent("""\
            Random number generator seed (-1 = arbitrary value). [%s]
            """ %('%(default)d')))

    g.add_argument('-pp', '--pc-permutations', type = int,
            metavar = int_mv, help = 'Number of permutations. [%s]'
            %('%(default)d'))

    g.add_argument('-pz', '--pc-zscore-thresh', type = float,
            metavar = float_mv, help = 'Z-score threshold. [%s]'
            %('%(default).2f'))

    g.add_argument('-pm', '--pc-max', type = int,
            metavar = int_mv, help = textwrap.dedent("""\
                Maximum number of PCs to test (0 = no maximum). [%s]
                """ %('%(default)d')))

    #for p in GOPCAConfig.param_defaults:
    #    assert p in dir(args)
    parser.set_defaults(**GOPCAConfig.param_defaults)

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
    config = None

    if args is None:
        # read arguments from the command line
        parser = get_argument_parser()

        # parse first with default options, in case "--help" is specified
        # ("--help" causes the program to exit at this point)
        args = parser.parse_args()

        # now remove the defaults and parse again
        # (removing the defaults is important so that we know which values
        # were specified by the user)
        no_defaults = dict([p, None] for p in GOPCAConfig.param_defaults)
        parser.set_defaults(**no_defaults)
        args = parser.parse_args()

    # reporting options
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # test if we can write to log_file?

    # configure root logger
    logger = util.get_logger(log_file = log_file, quiet = quiet,
            verbose = verbose)

    # generate configuration
    config = None

    if args.config_file is not None:
        # read parameter values from config file
        config = GOPCAConfig.read_config(args.config_file)

    if config is None:
        # start with default configuration
        config = GOPCAConfig()

    # overwrite parameters specified on the command line
    for p in GOPCAConfig.param_names:
        v = getattr(args, p)
        if v is not None:
            logger.debug('Parameter "%s" specified on command line!', p)
            config.set_param(p, v)

    # check if required parameters were specified
    passed = True
    if config.expression_file is None:
        logger.error('No expression file specified!')
        passed = False
    if config.gene_set_file is None:
        logger.error('No gene set file specified!')
        passed = False
    if config.output_file is None:
        logger.error('No output file specified!')
        passed = False
    if not passed:
        logger.error('Not all required parameters were specified.')
        return 1

    M = GOPCA(config)
    R = M.run()

    if R is None:
        logger.error('GO-PCA run failed!')
        return 1

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
