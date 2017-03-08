#!/usr/bin/env python

# Copyright (c) 2015-2017 Florian Wagner
#
# This file is part of GO-PCA.
#
# GO-PCA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Script to extract the GO-PCA signature matrix as a tab-delimited text file.
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
# import argparse
import logging

# import numpy as np
# from scipy.spatial.distance import pdist, squareform
# from scipy.cluster.hierarchy import linkage, dendrogram

from genometools import expression
from genometools import misc
from genometools.expression import ExpMatrix
from genometools.expression import cluster

from .. import util
from . import arguments


def get_argument_parser():

    desc = 'Extract the GO-PCA signatures matrix as a tab-delimited text file.'
    parser = arguments.get_argument_parser(desc=desc)

    arguments.add_io_args(parser)
    #arguments.add_signature_args(parser)
    #arguments.add_sample_args(parser)
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

    #sig_max_len = args.sig_max_len
    #sig_reverse_order = args.sig_reverse_order

    #sample_cluster_metric = args.sample_cluster_metric
    #no_sample_clustering = args.no_sample_clustering

    # configure root logger
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose
    logger = misc.get_logger(log_file=log_file, quiet=quiet,
                             verbose=verbose)

    result = util.read_gopca_result(gopca_file)
    
    sig_matrix = util.read_gopca_result(gopca_file)

    sig_labels = [sig.get_label(include_id=False)
                  for sig in sig_matrix.signatures]

    matrix = ExpMatrix(genes=sig_labels, samples=sig_matrix.samples,
                       X=sig_matrix.X)
    matrix.index.name = 'Signatures'
    #signatures = result.signatures
    #sig_labels = [sig.get_label(max_name_length=sig_max_len, include_id=False)
    #              for sig in signatures]
    #samples = list(result.samples)

    # generate expression matrix
    #E = ExpMatrix(genes=sig_labels, samples=samples, X=sig_matrix.X)

    # clustering of signatures (rows)
    #E, _ = cluster.cluster_genes(E, reverse=sig_reverse_order)

    exp_logger = logging.getLogger(expression.__name__)
    exp_logger.setLevel(logging.WARNING)
    matrix.write_tsv(output_file)
    exp_logger.setLevel(logging.NOTSET)
    logger.info('Wrote %d x %d signature matrix to "%s".',
                matrix.p, matrix.n, output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
