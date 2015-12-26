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

"""This script plots the GO-PCA signature matrix.

Example
-------

::

    $ gopca_plot_signature_matrix.py -g [gopca_output_file] -o [output_file]

"""

import sys
import os
import argparse
import cPickle as pickle
import logging

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

from genometools import misc

import gopca
from gopca import util
from gopca import cli
from gopca.plotting import cli as plot_cli

def get_argument_parser():

    prog = 'gopca_plot_signature_matrix.py'
    description = 'Plot the GO-PCA signature matrix.'
    parser = cli.get_argument_parser(prog, description)

    cli.add_io_args(parser)

    cli.add_reporting_args(parser)

    plot_cli.add_fig_args(parser)
    plot_cli.add_heatmap_args(parser)
    cli.add_signature_args(parser)
    cli.add_sample_args(parser)

    return parser

def main(args=None):

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    # input and output files
    gopca_file = args.gopca_file
    output_file = args.output_file

    # matplotlib backend
    mpl_backend = args.fig_mpl_backend

    # figure size
    fig_size = args.fig_size
    fig_res = args.fig_resolution

    # figure text
    use_tex = args.fig_use_tex
    font_size = args.fig_font_size
    font_family = args.fig_font_family

    # figure heatmap
    cmap = args.colormap
    vmin = args.val_coolest
    vmax = args.val_hottest

    # figure colorbar
    cbar_orient = args.cbar_orient
    cbar_anchor = args.cbar_anchor
    cbar_scale = args.cbar_scale
    cbar_pad = args.cbar_pad

    # samples
    sample_no_clustering = args.sample_no_clustering
    sample_cluster_metric = args.sample_cluster_metric

    # signatures
    sig_reverse_order = args.sig_reverse_order
    sig_max_len = args.sig_max_len
    sig_filter_corr = args.sig_filter_corr

    # reporting parameters
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # configure root logger
    logger = util.get_logger(log_file = log_file, quiet = quiet,
            verbose = verbose)

    # read GO-PCA output
    G = util.read_gopca_result(gopca_file)

    signatures = G.signatures
    S = G.S

    if sig_filter_corr < 1.0:
        q_before = G.q
        signatures, S = util.filter_signatures(signatures, S, sig_filter_corr)
        q = len(signatures)
        logger.info('Filtered %d / %d signatures.', q_before - q, q_before)
        
        #signature = sorted(G.signatures, 

    # generate labels
    labels = [sig.get_label(include_id=False,max_name_length=sig_max_len) for sig in signatures]
    samples = G.samples

    # clustering of rows (signatures)
    order_rows = util.cluster_signatures(S, reverse=sig_reverse_order)
    S = S[order_rows,:]
    labels = [labels[idx] for idx in order_rows]

    if not sample_no_clustering:
        # clustering of columns (samples)
        logger.info('Clustering of samples...')
        order_cols = util.cluster_samples(S, metric = sample_cluster_metric)
        S = S[:,order_cols]
        samples = [samples[j] for j in order_cols]

    # plotting
    logger.info('Plotting...')

    import matplotlib as mpl
    if mpl_backend is not None:
        mpl.use(mpl_backend)
    import matplotlib.pyplot as plt
    from matplotlib import rc

    #if plot_in_notebook:
    #   from IPython import get_ipython
    #   ipython = get_ipython()
    #   ipython.magic('matplotlib inline')

    if use_tex:
        rc('text', usetex=True)
    rc('font', family = font_family, size=font_size)
    rc('figure', figsize = (fig_size[0], fig_size[1]))
    rc('savefig', dpi = fig_res)

    # plotting
    plt.imshow(S, interpolation='none', aspect='auto',
            vmin=vmin, vmax=vmax, cmap=cmap)

    plt.xticks(())

    minint = int(vmin)
    maxint = int(vmax)
    cbticks = np.arange(minint, maxint+0.01, 1.0)
    cb = plt.colorbar(orientation = cbar_orient, shrink = cbar_scale,
            pad = cbar_pad, ticks=cbticks, use_gridspec=False,
            anchor=cbar_anchor)
    cb.ax.tick_params(labelsize = 'small')
    cb.set_label('Standardized Expression', size='small')

    q,n = S.shape
    plt.yticks(np.arange(q),labels,size='x-small')
    plt.xlabel('Samples (n=%d)' %(n))
    plt.ylabel('Signatures')

    logger.info('Saving to file...')
    plt.savefig(output_file,bbox_inches='tight')
    logger.info('Done!')

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
