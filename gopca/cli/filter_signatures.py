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

"""Script to filter redundant GO-PCA signatures.
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
# import cPickle as pickle
import textwrap

# import numpy as np

# from genometools import misc
from gopca import util
from gopca.cli import arguments
from gopca import GOPCASignatureMatrix


def get_argument_parser():

    desc = 'Filter redundant GO-PCA signatures.'
    parser = arguments.get_argument_parser(desc=desc)

    g = parser.add_argument_group('Filtering options')

    arguments.add_io_args(parser)
    float_mv = arguments.float_mv

    g.add_argument(
        '-r', '--corr-thresh', type=float, required=True, metavar=float_mv,
        help=textwrap.dedent("""\
            Correlation threshold for filtering signatures
            (1.0 = off)."""))

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
    output_file = args.output_file

    corr_thresh = args.corr_thresh

    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # configure root logger
    logger = util.get_logger(log_file=log_file, quiet=quiet,
                             verbose=verbose)

    result = util.read_gopca_result(gopca_file)

    signatures = result.signatures
    S = result.S

    if corr_thresh < 1.0:
        q_before = result.q
        signatures, S = util.filter_signatures(signatures, S, corr_thresh)
        q = len(signatures)
        logger.info('Filtered %d / %d signatures.', q_before-q, q_before)

    filtered = GOPCASignatureMatrix(result.config, result.genes, result.samples,
                                    result.W, result.Y, signatures, S)
    filtered.write_pickle(output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
