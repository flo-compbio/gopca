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

"""Script to print information about GO-PCA run and output data.
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
import argparse
import cPickle as pickle

import numpy as np

from genometools import misc
from genometools.expression import ExpMatrix
from gopca import util
from gopca import cli
from gopca import GOPCA, GOPCARun, GOPCAResult

def get_argument_parser():

    desc = 'Print information about GO-PCA run and output data.'
    parser = cli.get_argument_parser(desc = desc)

    str_type = cli.str_type
    file_mv = cli.file_mv

    g = parser.add_argument_group('Input file (required)')

    g.add_argument('-g', '--gopca-file', required = True,
            metavar = file_mv, type = str_type,
            help = 'A GO-PCA run or result pickle.')

    g.add_argument('-u', '--print-user-config', action = 'store_true',
            help = 'Print user-provided GO-PCA config data of the run.')

    g.add_argument('-s', '--print-signatures', action = 'store_true',
            help = 'Print signatures of the GO-PCA result.')

    cli.add_reporting_args(parser)

    return parser

def main(args=None):

    vinfo = sys.version_info
    if not (vinfo >= (2, 7)):
        raise SystemError('Python interpreter version >= 2.7 required, '
                          'found %d.%d instead.' %(vinfo.major, vinfo.minor))

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    gopca_file = args.gopca_file

    print_user_config = args.print_user_config
    print_signatures = args.print_signatures

    # configure root logger
    log_stream = sys.stderr
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    logger = util.get_logger(log_stream = log_stream, log_file = log_file,
            quiet = quiet, verbose = verbose)

    R = None
    G = None
    with open(gopca_file, 'rb') as fh:
        R = pickle.load(fh)

    if isinstance(R, GOPCARun):
        G = R.result
    elif isinstance(R, GOPCAResult):
        G = R
        R = None
    else:
        logger.error('The pickle contains neither GO-PCA run nor output data.')
        return 1

    if R is not None:
        print 'GO-PCA Run'
        print '----------'
        print '- GO-PCA version: %s' %(R.version)
        print '- Timestamp: %s' %(R.timestamp)
        print '- Exec. time: %.1f s' %(R.exec_time)

        if print_user_config:
            print '- User-provided config data:'
            for s in R.user_config.get_param_strings():
                print '    %s' %(s)
        print

    print 'GO-PCA Result'
    print '-------------'
    print '- Config MD5 hash: %s' %(G.config.get_hash())
    print '- Result MD5 hash: %s' %(G.get_hash())
    print '- Expression data: %d genes, %d samples' %(G.p, G.n)
    print '- Number of PCs tested: %d' %(G.D)
    print '- Number of signatures generated: %d' %(G.q)
    print '- Config data:'
    for s in G.config.get_param_strings():
        print '    %s' %(s)

    if print_signatures:
        GOPCA.print_signatures(G.signatures)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
