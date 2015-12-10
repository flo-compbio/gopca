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

"""This script plots a particular GO-PCA signature and its genes.

Example
-------

::

    $ gopca_plot_signature.py -g [gopca_output_file] -n [signature_name] -o [output_file]

"""

import sys
import os
import argparse
import cPickle as pickle

import numpy as np

from genometools import misc

import gopca
from gopca import util
from gopca import params
from gopca.plotting import params as plot_params

def get_argument_parser():

    prog = 'plot_gopca_signature.py'
    description = 'Plot a GO-PCA signature.'
    parser = params.get_argument_parser(prog, description)

    g = parser.add_argument_group('Required parameters')

    g.add_argument('-g', '--gopca-file', required = True,
            metavar = params.file_mv,
            help = 'The GO-PCA output file.')

    g.add_argument('-n', '--sig-name', required = True,
            metavar = params.name_mv,
            help = 'The name of the signature.')

    g.add_argument('-o', '--output-file', required = True,
            metavar = params.file_mv,
            help = 'The output file.')

    g = parser.add_argument_group('Layout')

    g.add_argument('-p', '--fig-title-pos', type = float, default = 0.95,
            metavar = params.float_mv,
            help = 'The position of the figure title.')

    g.add_argument('--fig-subgrid-ratio', type = int, default = 10,
            metavar = params.int_mv,
            help = 'The size ratio between signature and heat map panels.')

    g.add_argument('-gs', '--gene-label-size', type = float, default = None,
            metavar = params.float_mv,
            help = 'The size of the gene labels (in pt).')

    g.add_argument('-gr', '--gene-reverse-order', action = 'store_true',
            help = 'Reverse the order of the genes.')

    g.add_argument('--hide-id', action = 'store_true',
            help = 'Do not show the ID of the GO term.')

    plot_params.add_fig_params(parser)
    plot_params.add_heatmap_params(parser)
    params.add_sample_params(parser)

    return parser


def read_args_from_cmdline():

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-g','--gopca-file',required=True)
    parser.add_argument('-s','--sig-name',required=True)
    parser.add_argument('-o','--output-file',required=True)


    parser.add_argument('-d','--figure-dimensions',type=float,help='in inches',nargs=2,default=[15,10])
    parser.add_argument('-r','--figure-resolution',type=int,help='in dpi',default=150)
    parser.add_argument('-f','--figure-font-size',type=int,help='in pt',default=20)
    parser.add_argument('-m','--figure-font-family',default='serif')
    parser.add_argument('-c','--figure-colormap',default='RdBu_r')
    parser.add_argument('-vn','--figure-vmin',type=float,default=-3.0)
    parser.add_argument('-vx','--figure-vmax',type=float,default=3.0)

    parser.add_argument('-co','--figure-colorbar-orientation',default='horizontal')
    parser.add_argument('-ca','--figure-colorbar-anchor',type=float,nargs=2,default=(0.96,1.0))
    parser.add_argument('-cs','--figure-colorbar-shrink',type=float,default=0.3)
    parser.add_argument('-cp','--figure-colorbar-pad',type=float,default=0.015)


    parser.add_argument('-t','--use-tex',action='store_true')
    parser.add_argument('-b','--matplotlib-backend',default=None)
    parser.add_argument('-i','--invert-gene-order',action='store_true')
    parser.add_argument('--disable-sample-clustering',action='store_true')
    parser.add_argument('--cluster-global',action='store_true')
    parser.add_argument('--sort-by-signature',action='store_true')

    return parser.parse_args()

def main(args=None):

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    # required arguments
    gopca_file = args.gopca_file
    sig_name = args.sig_name
    output_file = args.output_file

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
    output = util.read_gopca_output(gopca_file)

    # find signature selected
    signatures = output.signatures
    term_ids = set([sig.term[0] for sig in signatures])
    sig = None
    if sig_name in term_ids:
        sig = [s for s in signatures if s.term[0] == sig_name]
        assert len(sig) == 1
        sig = sig[0]
    else:
        sig = [s for s in signatures
                if s.term[3].lower().startswith(sig_name.lower())]
        if len(sig) == 0:
            logger.error('Error: signature name "%s" not found.', sig_name)
            return 1
        elif len(sig) > 1:
            logger.error('Error: signature name not unique, matched: %s',
                    ', '.join([s.term[3] for s in sig]))
            return 1
        sig = sig[0]

    # get signature gene expression matrix and cluster rows
    sig_genes = sig.genes
    E = sig.E
    logger.debug('Expression matrix shape', str(E.shape))
    order_rows = util.cluster_rows(E, metric='correlation', method='average')
    if gene_reverse_order:
        order_rows = order_rows[::-1]
    sig_genes = [sig_genes[i] for i in order_rows]
    E = E[order_rows,:]

    # standardize gene expression matrix
    E_std = util.get_standardized_matrix(E)
    # calculate signature label and expression
    include_id = not(hide_id)
    sig_label = sig.get_label(include_id = include_id)
    sig_expr = np.mean(E_std, axis=0)

    # get sample order from signature matrix
    #if sort_by_signature:
    #    a = np.argsort(sig_expr)
    #    a = a[::-1]
    #    sig_expr = sig_expr[a]
    #    E_std = E_std[:,a]
    if (not sample_no_clustering):
        S = output.S
        order_cols = util.cluster_samples(S, metric = sample_cluster_metric)
        sig_expr = sig_expr[order_cols]
        E_std = E_std[:,order_cols]

    # plotting
    logger.info('Plotting...')

    import matplotlib as mpl
    if mpl_backend is not None:
        mpl.use(mpl_backend)
    import matplotlib.pyplot as plt
    from matplotlib import rc

    if use_tex: rc('text', usetex=True)
    rc('font', family = font_family, size = font_size)
    rc('figure', figsize = fig_size)
    rc('savefig', dpi = fig_res)

    # subgrid layout
    ax = plt.subplot2grid((subgrid_ratio, 1), (0, 0))
    plt.sca(ax)

    plt.imshow(np.atleast_2d(sig_expr), aspect = 'auto',\
            interpolation = 'none', vmin = vmin, vmax = vmax, cmap = cmap)
    plt.xticks(())
    plt.yticks([0], ['Signature'] , size='small')

    ax = plt.subplot2grid((subgrid_ratio, 1), (1, 0),
            rowspan = subgrid_ratio - 1)
    plt.sca(ax)
    plt.imshow(E_std, interpolation='none', aspect='auto',
            vmin = vmin, vmax = vmax, cmap = cmap)
    plt.yticks(np.arange(len(sig_genes)), sig_genes , size = gene_label_size)
    plt.xlabel('Samples')
    plt.ylabel('Genes')
    plt.xticks(())

    minint = int(vmin)
    maxint = int(vmax)
    cbticks = np.arange(minint, maxint + 0.01, 1.0)
    cb = plt.colorbar(orientation = cbar_orient, shrink = cbar_scale,
            pad = cbar_pad, ticks = cbticks, use_gridspec = False,
            anchor = cbar_anchor)
    cb.ax.tick_params(labelsize = 'small')
    cb.set_label('Standardized Expression',size = 'small')
    plt.suptitle(sig_label, va = 'top', y = title_pos)

    #plt.tight_layout()

    logger.info('Saving to file...')
    plt.savefig(output_file, bbox_inches='tight')
    logger.info('Done!')

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
