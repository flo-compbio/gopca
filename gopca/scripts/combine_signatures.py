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

"""Script to combine the signatures from multiple GO-PCA analyses.

Assumes that the sample names are identical in all analyses.
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
import cPickle as pickle
import textwrap

import numpy as np

from genometools import misc
from gopca import util
from gopca import cli
#from gopca import GOPCAResult

def get_argument_parser():

    desc ='Combine the signatures from multiple GO-PCA analyses.'
    parser = cli.get_argument_parser(desc = desc)

    str_type = cli.str_type
    file_mv = cli.file_mv

    g = parser.add_argument_group('Input and output files')

    g.add_argument('-g', '--gopca-files', nargs = '+', required = True,
            type = str_type, metavar = file_mv,
            help = 'List of GO-PCA output files with signatures to be merged.')

    g.add_argument('-o', '--output-file', required = True,
            type = str_type, metavar = file_mv,
            help = 'The output file.')

    cli.add_reporting_args(parser)

    return parser

def main(args = None):

    vinfo = sys.version_info
    if not (vinfo >= (2, 7)):
        raise SystemError('Python interpreter version >= 2.7 required, '
                          'found %d.%d instead.' %(vinfo.major, vinfo.minor))

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

    # read first result
    G = util.read_gopca_result(gopca_files[0])

    # set all GOPCAResult attributes to None,
    # except ``samples``, ``S`` and ``signatures``
    G.config = None
    G.genes = None
    G.W = None
    G.Y = None

    # append other outputs
    for i in range(1, len(gopca_files)):
        G_other = util.read_gopca_result(gopca_files[i])
        assert tuple(G.samples) == tuple(G_other.samples)

        G.signatures = G.signatures + G_other.signatures
        G.S = np.vstack([G.S, G_other.S])

    G.write_pickle(output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
