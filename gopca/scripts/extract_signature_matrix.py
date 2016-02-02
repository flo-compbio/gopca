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

"""Script to extract the GO-PCA signature matrix as a tab-delimited text file.
"""

import sys
import argparse
import logging

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

from genometools import expression
from genometools import misc
from genometools.expression import ExpMatrix
from gopca import util
from gopca import cli

def get_argument_parser():

    desc = 'Extract the GO-PCA signatures matrix as a tab-delimited text file.'
    parser = cli.get_argument_parser(desc = desc)

    cli.add_io_args(parser)
    cli.add_signature_args(parser)
    cli.add_sample_args(parser)
    cli.add_reporting_args(parser)

    return parser

def main(args = None):

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    gopca_file = args.gopca_file
    output_file = args.output_file

    sig_max_len = args.sig_max_len
    sig_reverse_order = args.sig_reverse_order

    sample_cluster_metric = args.sample_cluster_metric
    sample_no_clustering = args.sample_no_clustering

    # configure root logger
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    logger = misc.get_logger(log_file = log_file, quiet = quiet,
            verbose = verbose)

    G = util.read_gopca_result(gopca_file)
    
    signatures = G.signatures
    labels = [sig.get_label(include_id=False) for sig in signatures]
    samples = list(G.samples)
    S = G.S

    # clustering of rows (signatures)
    order_rows = util.cluster_signatures(S, reverse = sig_reverse_order)
    S = S[order_rows,:]
    labels = [labels[idx] for idx in order_rows]

    if not sample_no_clustering:
        # clustering of columns (samples)
        #distxy = squareform(pdist(S.T, metric='euclidean'))
        distxy = squareform(pdist(S.T, metric = sample_cluster_metric))
        R = dendrogram(linkage(distxy, method='average'), no_plot=True)
        order_cols = np.int64([int(l) for l in R['ivl']])
        S = S[:,order_cols]
        samples = [samples[j] for j in order_cols]

    exp_logger = logging.getLogger(expression.__name__)
    exp = ExpMatrix(labels, samples, S)
    exp_logger.setLevel(logging.WARNING)
    exp.write_tsv(output_file)
    exp_logger.setLevel(logging.NOTSET)
    logger.info('Wrote matrix with %d signatures and %d samples to "%s".',
            len(signatures), len(samples), output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
