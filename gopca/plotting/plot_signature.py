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

import sys
import os
import argparse
import cPickle as pickle

import numpy as np

from gopca import common
from genometools import misc

def read_args_from_cmdline():

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-g','--gopca-file',required=True)
    parser.add_argument('-s','--signature-name',required=True)
    parser.add_argument('-o','--output-file',required=True)

    parser.add_argument('-l','--gene-label-size',default='x-small')

    parser.add_argument('-d','--figure-dimensions',type=float,help='in inches',nargs=2,default=[15,10])
    parser.add_argument('-r','--figure-resolution',type=int,help='in dpi',default=150)
    parser.add_argument('-f','--figure-font-size',type=int,help='in pt',default=20)
    parser.add_argument('-m','--figure-font-family',default='serif')
    parser.add_argument('-c','--figure-colormap',default='RdBu_r')
    parser.add_argument('-vn','--figure-vmin',type=float,default=-3.0)
    parser.add_argument('-vx','--figure-vmax',type=float,default=3.0)
    parser.add_argument('-p','--figure-title-pos',type=float,default=0.95)

    parser.add_argument('-co','--figure-colorbar-orientation',default='horizontal')
    parser.add_argument('-ca','--figure-colorbar-anchor',type=float,nargs=2,default=(0.96,1.0))
    parser.add_argument('-cs','--figure-colorbar-shrink',type=float,default=0.3)
    parser.add_argument('-cp','--figure-colorbar-pad',type=float,default=0.015)

    parser.add_argument('--exclude-id',action='store_true')

    parser.add_argument('--figure-subgrid-ratio',type=int,default=10)

    parser.add_argument('-t','--use-tex',action='store_true')
    parser.add_argument('-b','--matplotlib-backend',default=None)
    parser.add_argument('-i','--invert-gene-order',action='store_true')
    parser.add_argument('--disable-sample-clustering',action='store_true')
    parser.add_argument('--cluster-global',action='store_true')
    parser.add_argument('--sort-by-signature',action='store_true')

    return parser.parse_args()

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def main(args=None):

    if args is None:
        args = read_args_from_cmdline()

    result_file = args.gopca_file
    sig_name = args.signature_name
    output_file = args.output_file

    gene_label_size = args.gene_label_size
    include_id = not(args.exclude_id)

    # figure size
    fig_dim = args.figure_dimensions
    fig_res = args.figure_resolution

    # figure text
    use_tex = args.use_tex
    fig_font_size = args.figure_font_size
    fig_font_family = args.figure_font_family

    # figure heatmap
    fig_vmin = args.figure_vmin
    fig_vmax = args.figure_vmax
    fig_cmap = args.figure_colormap

    # figure colorbar
    fig_cbar_orient = args.figure_colorbar_orientation
    fig_cbar_anchor = args.figure_colorbar_anchor
    fig_cbar_shrink = args.figure_colorbar_shrink
    fig_cbar_pad = args.figure_colorbar_pad

    # specific
    fig_title_pos = args.figure_title_pos
    fig_subgrid_ratio = args.figure_subgrid_ratio
    disable_sample_clustering = args.disable_sample_clustering
    cluster_global = args.cluster_global
    sort_by_signature = args.sort_by_signature

    mpl_backend = args.matplotlib_backend
    invert_gene_order = args.invert_gene_order

    # read GO-PCA result
    result = None
    with open(result_file,'rb') as fh:
        result = pickle.load(fh)

    # find signature selected
    signatures = result.signatures
    term_ids = set([sig.term[0] for sig in signatures])
    sig = None
    if sig_name in term_ids:
        sig = [s for s in signatures if s.term[0] == sig_name]
        assert len(sig) == 1
        sig = sig[0]
    else:
        sig_name = sig_name.lower()
        sig = [s for s in signatures if s.term[3].lower().startswith(sig_name)]
        if len(sig) == 0:
            print >> sys.stderr, 'Error: signature name not found.'
            return 1
        elif len(sig) > 1:
            print >> sys.stderr, 'Error: signature name not unique, matched: %s' %(', '.join([s.term[3] for s in sig]))
            return 1
        sig = sig[0]

    # get signature gene expression matrix and cluster rows
    sig_genes = sig.genes
    E = sig.E
    print E.shape
    order_rows = common.cluster_rows(E,metric='correlation',method='average')
    if invert_gene_order:
        order_rows = order_rows[::-1]
    sig_genes = [sig_genes[i] for i in order_rows]
    E = E[order_rows,:]

    # standardize gene expression matrix
    E_std = common.get_standardized_matrix(E)
    # calculate signature label and expression
    sig_label = sig.get_label(include_id=include_id)
    sig_expr = np.mean(E_std,axis=0)


    # get sample order from signature matrix
    if sort_by_signature:
        a = np.argsort(sig_expr)
        a = a[::-1]
        sig_expr = sig_expr[a]
        E_std = E_std[:,a]
    elif not disable_sample_clustering:
        if cluster_global:
            S = result.S
            order_cols = common.cluster_samples(S)
        else:
            order_cols = common.cluster_samples(E)
        sig_expr = sig_expr[order_cols]
        E_std = E_std[:,order_cols]

    # plotting
    import matplotlib as mpl
    if mpl_backend is not None:
        mpl.use(mpl_backend)
    import matplotlib.pyplot as plt
    from matplotlib import rc

    if use_tex: rc('text',usetex=True)
    rc('font',family=fig_font_family,size=fig_font_size)
    rc('figure',figsize=(fig_dim[0],fig_dim[1]))
    rc('savefig',dpi=fig_res)

    # subgrid layout
    ax = plt.subplot2grid((fig_subgrid_ratio,1),(0,0))
    plt.sca(ax)

    plt.imshow(np.atleast_2d(sig_expr),aspect='auto',interpolation='none',vmin=fig_vmin,vmax=fig_vmax,cmap=fig_cmap)
    plt.xticks(())
    plt.yticks([0],['Signature'],size='small')

    ax = plt.subplot2grid((fig_subgrid_ratio,1),(1,0),rowspan=fig_subgrid_ratio-1)
    plt.sca(ax)
    plt.imshow(E_std,interpolation='none',aspect='auto',vmin=fig_vmin,vmax=fig_vmax,cmap=fig_cmap)
    plt.yticks(np.arange(len(sig_genes)),sig_genes,size=gene_label_size)
    plt.xlabel('Samples')
    plt.ylabel('Genes')
    plt.xticks(())

    minint = int(fig_vmin)
    maxint = int(fig_vmax)
    cbticks = np.arange(minint,maxint+0.01,1.0)
    cb = plt.colorbar(orientation=fig_cbar_orient,shrink=fig_cbar_shrink,pad=fig_cbar_pad,ticks=cbticks,use_gridspec=False,anchor=fig_cbar_anchor)
    cb.ax.tick_params(labelsize='small')
    cb.set_label('Standardized Expression',size='small')
    plt.suptitle(sig_label,va='top',y=fig_title_pos)

    #plt.tight_layout()

    print 'Saving to file...', ; sys.stdout.flush()
    plt.savefig(output_file,bbox_inches='tight')
    print 'done!'; sys.stdout.flush()

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
