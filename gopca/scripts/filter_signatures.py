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

"""Script to filter redundant GO-PCA signatures.
"""

import sys
import cPickle as pickle
import textwrap

import numpy as np

from genometools import misc
from gopca import util
from gopca import cli
from gopca import GOPCAResult

def get_argument_parser():

    desc ='Filter redundant GO-PCA signatures.'
    parser = cli.get_argument_parser(desc = desc)

    g = parser.add_argument_group('Filtering options')

    cli.add_io_args(parser)

    g.add_argument('-r', '--corr-thresh', type = float, required = True,
            metavar = cli.float_mv,
            help = textwrap.dedent("""\
                Correlation threshold for filtering signatures
                (1.0 = off)."""))

    cli.add_reporting_args(parser)

    return parser

def main(args = None):

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
    logger = util.get_logger(log_file = log_file, quiet = quiet,
            verbose = verbose)

    G = util.read_gopca_result(gopca_file)

    signatures = G.signatures
    S = G.S

    if corr_thresh < 1.0:
        q_before = G.q
        signatures, S = util.filter_signatures(signatures, S, corr_thresh)
        q = len(signatures)
        logger.info('Filtered %d / %d signatures.', q_before - q, q_before)

    F = GOPCAResult(G.config, G.genes, G.samples, G.W, G.Y, signatures, S)
    F.write_pickle(output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
