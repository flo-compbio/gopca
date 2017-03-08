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
# import argparse

import six
# import numpy as np

# from genometools import misc
# from genometools.expression import ExpMatrix
from .. import util
from .. import GOPCA, GOPCARun, GOPCASignatureMatrix
from . import arguments

if six.PY2:
    import cPickle as pickle
else:
    import pickle


def get_argument_parser():

    desc = 'Print information about GO-PCA run and output data.'
    parser = arguments.get_argument_parser(desc=desc)

    file_mv = arguments.file_mv

    g = parser.add_argument_group('Input file (required)')

    g.add_argument(
        '-g', '--gopca-file', type=str, required=True, metavar=file_mv,
        help='A GO-PCA run or result pickle.')

    g.add_argument(
        '-u', '--print-user-config', action='store_true',
        help='Print user-provided GO-PCA config data of the run.')

    g.add_argument(
        '-s', '--print-signatures', action='store_true',
        help='Print signatures of the GO-PCA result.')

    arguments.add_reporting_args(parser)

    return parser


def main(args=None):

    vinfo = sys.version_info
    if not (vinfo >= (2, 7)):
        raise SystemError('Python interpreter version >= 2.7 required, '
                          'found %d.%d instead.' % (vinfo.major, vinfo.minor))

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

    logger = util.get_logger(log_stream=log_stream, log_file=log_file,
                             quiet=quiet, verbose=verbose)

    with open(gopca_file, 'rb') as fh:
        run = pickle.load(fh)

    if isinstance(run, GOPCARun):
        result = run.sig_matrix
    elif isinstance(run, GOPCASignatureMatrix):
        result = run
        run = None
    else:
        logger.error('The pickle contains neither GO-PCA run nor output data.')
        return 1

    if run is not None:
        print('GO-PCA Run')
        print('----------')
        print('- GO-PCA version: %s' % run.gopca_version)
        print('- Timestamp: %s' % run.timestamp)
        print('- Exec. time: %.1f s' % run.exec_time)

        print('- Expression data MD5 hash: %s' % run.expression_hash)
        print('- Number of configurations: %d' % len(run.config_hashes))
        for i, ch in enumerate(run.config_hashes):
            print('\tConfiguration %d MD5 hash: %s' % (i+1, ch))

        #if print_user_config:
        #    print('- User-provided config data:')
        #    for s in run.user_config.get_param_strings():
        #        print('    %s' % s)
        print()

    print('GO-PCA Signature Matrix')
    print('-----------------------')
    print('- MD5 hash: %s' % result.hash)
    #print('- Config MD5 hash: %s' % result.config.get_hash())
    #print('- Result MD5 hash: %s' % result.get_hash())
    #print('- Expression data: %d genes, %d samples' % (result.p, result.n))
    #print('- Number of PCs tested: %d' % result.D)
    print('- Number of signatures generated: %d' % len(result.index))
    #print('- Config data:')
    #for s in result.config.get_param_strings():
    #    print('    %s' % s)

    #if print_signatures:
    #    GOPCA.print_signatures(result.signatures)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
