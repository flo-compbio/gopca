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

"""Script to combine multiple GO-PCA outputs.
"""

import sys
import cPickle as pickle
import textwrap

import numpy as np

from genometools import misc
from gopca import util
from gopca import cli
from gopca import GOPCAOutput

def get_argument_parser():

    desc ='Combine multiple GO-PCA outputs.'
    parser = cli.get_argument_parser(description = desc)

    str_type = cli.str_type
    file_mv = cli.file_mv

    g = parser.add_argument_group('Input and output files')

    g.add_argument('-g', '--gopca-files', nargs = '+', required = True,
            type = str_type, metavar = file_mv)

    g.add_argument('-o', '--output-file', required = True,
            type = str_type, metavar = file_mv)

    cli.add_reporting_args(parser)

    return parser

def main(args = None):

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    gopca_files = args.gopca_files
    output_file = args.output_file

    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # configure root logger
    logger = util.get_logger(log_file = log_file, quiet = quiet,
            verbose = verbose)

    # read first output
    G = util.read_gopca_output(gopca_files[0])

    config = G.config
    signatures = G.signatures
    S = G.S
    genes = G.genes
    samples = G.samples
    W = G.W
    Y = G.Y

    # append other outputs
    for i in range(1, len(gopca_files)):
        G = util.read_gopca_output(gopca_files[i])
        assert tuple(G.genes) == tuple(genes)
        assert tuple(G.samples) == tuple(samples)

        signatures = signatures + G.signatures
        S = np.vstack([S,G.S])

    G = GOPCAOutput(config, genes, samples, W, Y, signatures, S)
    G.write_pickle(output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
