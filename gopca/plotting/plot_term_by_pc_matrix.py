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

"""This script plots the GO-PCA term-by-PC matrix.

The matrix shows associations between GO terms which were used to generate
signatures and principal components (PCs).

This visualization was originally proposed by Dr. Meromit Singer.

Example
-------

::

    $ gopca_plot_term_by_pc_matrix.py -g [gopca_output_file] -o [output_file]

"""


import sys
import os
import argparse
import csv
import cPickle as pickle
from math import ceil
from collections import OrderedDict,Counter

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

import matplotlib as mpl
#from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc

import xlmhg

from genometools import misc

import gopca
from gopca import util
from gopca import cli
from gopca.plotting import cli as plot_cli

def get_argument_parser():

    desc = 'Plot the GO-PCA term-by-PC matrix.'
    parser = cli.get_argument_parser(desc = desc)

    cli.add_io_args(parser)

    g = parser.add_argument_group('Optional parameters')

    g.add_argument('--max-pval', type = float, default = 1e-4)
    g.add_argument('--dotcolor', default = 'yellow')
    g.add_argument('--dotsize', type = float, default = 50)

    plot_cli.add_fig_args(parser)
    parser.set_defaults(fig_font_size = 16)

    plot_cli.add_heatmap_args(parser)
    parser.set_defaults(val_coolest = 0.0)
    parser.set_defaults(val_hottest = 10.0)
    parser.set_defaults(cbar_ticks = [4,6,8,10])
    parser.set_defaults(colormap = 'BuPu')
    parser.set_defaults(cbar_anchor = [0.6,1.0])
    parser.set_defaults(cbar_scale = 0.2)

    cli.add_go_term_args(parser)
    cli.add_sample_args(parser)

    return parser

def main(args=None):

    if args is None:
        # read command line arguments
        parser = get_argument_parser()
        args = parser.parse_args()

    # input and output files
    gopca_file = args.gopca_file
    output_file = args.output_file

    # options
    max_pval = args.max_pval
    dotcolor = args.dotcolor
    dotsize = args.dotsize

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
    cbar_ticks = args.cbar_ticks

    # GO terms
    term_reverse_order = args.term_reverse_order
    term_max_len = args.term_max_len

    # configure root logger
    #logger = misc.configure_logger('')

    # read GO-PCA output
    G = util.read_gopca_output(gopca_file)

    genes = G.genes
    signatures = G.signatures
    W = G.W # loading matrix
    S = G.S # signature matrix
    #input_ = G.input
    config = G.config
    X_frac = config.mHG_X_frac
    X_min = config.mHG_X_min
    L = config.mHG_L

    # order genes alphabetically
    a = np.lexsort([genes])
    p = len(genes)
    assert np.all(a == np.arange(p))
    genes = [genes[i] for i in a]
    W = W[a,:]

    # determine max K
    K_max = max(sig.K for sig in signatures)

    # generate term labels
    labels = ['%s (%d)' %(sig.get_label(max_name_length = term_max_len,
            include_id = False, include_stats = False), sig.K)
            for sig in signatures]

    # order terms using hierarchical clustering
    order = util.cluster_signatures(S)
    if term_reverse_order:
        order = order[::-1]
    S = S[order,:]
    labels = [labels[i] for i in order]
    signatures = [signatures[i] for i in order]

    # test association
    p,n_comps = W.shape
    q = len(signatures)
    A = np.zeros((q, 2*n_comps), dtype=np.float64)
    matrix = np.empty((K_max+1, p+1), dtype=np.longdouble)

    for pc in range(n_comps):
        a_asc = np.argsort(W[:,pc])
        a_dsc = a_asc[::-1]
        for i,sig in enumerate(signatures):
            tg = sig.genes # != sig.enr.genes
            K = len(tg)
            X = max(X_min,int(ceil(X_frac*K)))
        
            v = np.zeros(p,dtype=np.uint8)
            for g in tg:
                idx = misc.bisect_index(genes,g)
                v[idx] = 1

            v_sorted = np.ascontiguousarray(v[a_dsc])
            threshold,_,pval = xlmhg.test(v_sorted,X,L,K,mat=matrix)
            A[i,pc*2] = -np.log10(pval)

            v_sorted = np.ascontiguousarray(v[a_asc])
            threshold,_,pval = xlmhg.test(v_sorted,X,L,K,mat=matrix)
            A[i,pc*2+1] = -np.log10(pval)

    # prepare for plotting
    rc('font',family = font_family, size = font_size)
    rc('figure', figsize = fig_size)
    rc('savefig', dpi = fig_res)

    if use_tex:
        rc('text', usetex=True)
        preamble = mpl.rcParams['text.latex.preamble']
        add = r'\usepackage{bm}'
        if add not in preamble:
            mpl.rcParams['text.latex.preamble'].append(add)

    # for each signature, mark PC that was originally used to generate it
    # plot this first, otherwise colormap gets messed up
    q = len(signatures)
    for i in range(q):
        j = (abs(signatures[i].pc)-1)*2
        if signatures[i].pc < 0:
            j += 1
        plt.scatter([j], [i], color = dotcolor, zorder = 100, s = dotsize,
                marker = 'x')

    # plot heatmap
    A[np.absolute(A) < -np.log10(max_pval)] = np.nan # hide insignificant associations
    plt.imshow(A, interpolation = 'none', vmin = vmin, vmax = vmax, cmap = cmap, zorder=20)

    # plot separation lines
    for pc in range(n_comps-1):
        plt.plot([pc*2+1.5,pc*2+1.5],[-0.5,q-0.5],color='gray',linewidth=1.0,zorder=50)
    plt.plot([-0.5,n_comps*2+1.5],[-0.5,-0.5],color='black',linewidth=2.0,zorder=50) # fix some z-order issues

    # configure axes
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('top')
    plt.gca().xaxis.set_label_position('top')
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
 
    #plt.gca().xaxis.tick_top()
    plt.xticks(np.arange(0,n_comps*2,2)+0.5,np.arange(n_comps)+1,size='small')
    plt.xlabel(r'Principal Component',labelpad=10,size='small')
    plt.xlim(-0.5,n_comps*2-0.5)
    plt.gca().tick_params(top='off')

    plt.yticks(np.arange(q), labels, size='x-small')
    plt.ylabel(r'GO Term',size='small')
    plt.ylim(q-0.5,-0.5)
    #plt.ylim(-0.5,q-0.5)
    plt.grid(which='both',axis='y',zorder=-20) # z-order is ignored here

    # plot colorbar
    cb = plt.colorbar(orientation = cbar_orient, shrink = cbar_scale, pad = cbar_pad, ticks = cbar_ticks, use_gridspec = False, anchor = cbar_anchor)
    cb.ax.tick_params(labelsize = 'small')

    if use_tex:
        cb.set_label(r"$\bm{-\log_{10}} \textrm{p-value}$",size='small')
    else:
        cb.set_label('-Log10 p-value',size='small')
    plt.savefig(output_file,bbox_inches='tight')

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
