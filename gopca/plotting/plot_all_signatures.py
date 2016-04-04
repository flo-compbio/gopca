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

"""This script plots all GO-PCA signatures to a PDF file (one sig. per page).

Example
-------

::

    $ gopca_plot_all_signatures.py -g [gopca_output_file] -o [output_file]

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
import os
import argparse
import cPickle as pickle

import numpy as np

from genometools import misc

import gopca
from gopca import util
from gopca import cli
from gopca.plotting import cli as plot_cli

def get_argument_parser():

    desc = 'Plot a GO-PCA signature.'
    parser = cli.get_argument_parser(desc = desc)

    g = cli.add_io_args(parser)

    g = cli.add_signature_args(parser)

    g.add_argument('--no-standardization', action = 'store_true',
            help = 'Do not standardize the signature expression values.')

    g = parser.add_argument_group('Layout options')

    g.add_argument('-p', '--fig-title-pos', type = float, default = 0.95,
            metavar = cli.float_mv,
            help = 'The position of the figure title.')

    g.add_argument('--fig-subgrid-ratio', type = int, default = 10,
            metavar = cli.int_mv,
            help = 'The size ratio between signature and heat map panels.')

    g.add_argument('-gs', '--gene-label-size', type = float, default = None,
            metavar = cli.float_mv,
            help = 'The size of the gene labels (in pt).')

    g.add_argument('-gr', '--gene-reverse-order', action = 'store_true',
            help = 'Reverse the order of the genes.')

    g.add_argument('--hide-id', action = 'store_true',
            help = 'Do not show the ID of the GO term.')

    plot_cli.add_fig_args(parser)
    plot_cli.add_heatmap_args(parser)
    cli.add_sample_args(parser)

    return parser

def main(args = None):

    vinfo = sys.version_info
    if not (vinfo >= (2, 7)):
        raise SystemError('Python interpreter version >= 2.7 required, '
                          'found %d.%d instead.' %(vinfo.major, vinfo.minor))

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    # input/output
    gopca_file = args.gopca_file
    output_file = args.output_file

    # signature
    sig_reverse_order = args.sig_reverse_order
    no_standardization = args.no_standardization

    # layout
    gene_label_size = args.gene_label_size
    if gene_label_size is None:
        gene_label_size = 'x-small'

    hide_id = args.hide_id
    title_pos = args.fig_title_pos
    subgrid_ratio = args.fig_subgrid_ratio
    gene_reverse_order = args.gene_reverse_order
    #cluster_global = args.cluster_global
    #sort_by_signature = args.sort_by_signature

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
    sample_font_size = args.sample_label_font_size
    if sample_font_size is None:
        sample_font_size = int(0.75*font_size)

    # figure colorbar
    cbar_orient = args.cbar_orient
    cbar_anchor = args.cbar_anchor
    cbar_scale = args.cbar_scale
    cbar_pad = args.cbar_pad

    # samples
    sample_no_clustering = args.sample_no_clustering
    sample_cluster_metric = args.sample_cluster_metric

    # configure root logger
    logger = misc.configure_logger('')

    # read GO-PCA output
    G = util.read_gopca_result(gopca_file)
    samples = G.samples
    signatures = G.signatures
    S = G.S

    # clustering of rows (signatures)
    order_rows = util.cluster_signatures(S, reverse = sig_reverse_order)
    signatures = [signatures[i] for i in order_rows]
    S = S[order_rows,:]
    #labels = [labels[idx] for idx in order_rows]

    import matplotlib as mpl
    #if mpl_backend is not None:
    #    mpl.use(mpl_backend)
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    from matplotlib import rc
    
    with PdfPages(output_file) as pdf:

        for k, sig in enumerate(signatures):

            sig_label = sig.get_label(include_id = True)
            logger.info('Plotting signature: %s', sig_label)

            # get signature gene expression matrix and cluster rows
            sig_genes = sig.genes
            X = sig.X
            logger.debug('Expression matrix shape', str(X.shape))
            order_rows = util.cluster_rows(X, metric = 'correlation',
                    method = 'average')
            if gene_reverse_order:
                order_rows = order_rows[::-1]
            sig_genes = [sig_genes[i] for i in order_rows]
            X = X[order_rows,:]

            if (not no_standardization):
                # standardize gene expression matrix
                X_std = util.get_standardized_matrix(X)
            else:
                # only center expression matrix
                X_std = util.get_centered_matrix(X)
            # calculate signature label and expression
            include_id = not(hide_id)
            sig_label = sig.get_label(include_id = include_id)
            sig_expr = np.mean(X_std, axis=0)

            # get sample order from signature matrix
            #if sort_by_signature:
            #    a = np.argsort(sig_expr)
            #    a = a[::-1]
            #    sig_expr = sig_expr[a]
            #    E_std = E_std[:,a]
            sample_labels = samples
            if (not sample_no_clustering):
                S = G.S
                order_cols = util.cluster_samples(S, metric = sample_cluster_metric)
                sig_expr = sig_expr[order_cols]
                X_std = X_std[:,order_cols]
                sample_labels = [sample_labels[j] for j in order_cols]

            # plotting
            logger.info('Plotting...')

            if use_tex: rc('text', usetex=True)
            rc('font', family = font_family, size = font_size)
            rc('figure', figsize = fig_size)
            rc('savefig', dpi = fig_res)

            # erase figure
            plt.clf()
            plt.cla()
            #plt.gcf().subplots_adjust(left=0.30)
            #plt.gcf().subplots_adjust(right=0.95)

            #plt.subplots_adjust(left=None, bottom=0, right=None, top=None,
            #    wspace=None, hspace=0)

            # subgrid layout
            ax = plt.subplot2grid((subgrid_ratio, 1), (0, 0))
            plt.sca(ax)

            #plt.gcf().subplots_adjust(bottom = 0.1)
            #print plt.subplots_adjust()

            plt.imshow(np.atleast_2d(sig_expr), aspect = 'auto',\
                    interpolation = 'none', vmin = vmin, vmax = vmax, cmap = cmap)
            plt.xticks(())
            plt.yticks([0], ['Signature'] , size='small')

            ax = plt.subplot2grid((subgrid_ratio, 1), (1, 0),
                    rowspan = subgrid_ratio - 1)
            #print 'Position:', ax.get_position()
            #print 'Window extent', ax.get_window_extent()
            plt.sca(ax)

            plt.imshow(X_std, interpolation='none', aspect='auto',
                    vmin = vmin, vmax = vmax, cmap = cmap)
            q, n = S.shape
            plt.xticks(())
            if args.show_sample_labels:
                plt.xticks(np.arange(n), sample_labels, size = sample_font_size, rotation = 30, ha = 'right')
            plt.yticks(np.arange(len(sig_genes)), sig_genes , size = gene_label_size)
            plt.xlabel('Samples')
            plt.ylabel('Genes')

            minint = int(vmin)
            maxint = int(vmax)
            cbticks = np.arange(minint, maxint + 0.01, 1.0)
            if args.show_sample_labels:
                cbar_pad += 0.1
            cb = plt.colorbar(orientation = cbar_orient, shrink = cbar_scale,
                    pad = cbar_pad, ticks = cbticks, use_gridspec = False,
                    anchor = cbar_anchor)
            cb.ax.tick_params(labelsize = 'small')

            if (not no_standardization):
                cb.set_label('Standardized Expression',size = 'small')
            else:
                cb.set_label('Centered Expression', size = 'small')
            plt.suptitle(sig_label, va = 'top', y = title_pos)

            #plt.tight_layout()

            pdf.savefig(bbox_inches = 'tight')
            #plt.savefig(output_file, bbox_inches='tight')
            #if k >= 4:
            #    break

    logger.info('Saving to file...')
    logger.info('Done!')
    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
