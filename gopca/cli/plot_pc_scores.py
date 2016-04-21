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

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
# import os
import argparse

import numpy as np

# from genometools import misc
from gopca import util
from gopca.cli import arguments


def get_argument_parser():

    desc = 'Plot PC scores with signature values overlaid.'
    parser = arguments.get_argument_parser(desc=desc)

    g = arguments.add_io_args(parser)

    str_mv = arguments.str_mv
    int_mv = arguments.int_mv

    g = parser.add_argument_group('Signature options')

    g.add_argument(
        '-n', '--sig-name', type=str, required=True, metavar=str_mv,
        help='The GO term ID or name of the signature.')

    g.add_argument(
        '-c', '--components', type=int, nargs=2, default=[1, 2],
        metavar=int_mv,
        help='The principal components to show.')

    arguments.add_fig_args(parser)

    arguments.add_heatmap_args(parser)

    g = arguments.add_reporting_args(parser)

    """
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
    """

    return parser


def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


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

    sig_name = args.sig_name
    comps = args.components

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

    # reporting parameters
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # configure root logger
    logger = util.get_logger(log_file=log_file, quiet=quiet,
                             verbose=verbose)

    # read GO-PCA result
    result = util.read_gopca_result(gopca_file)

    # find signature selected
    term_ids = set([sig.term[0] for sig in result.signatures])
    # sig = None
    # idx = None
    if sig_name in term_ids:
        idx = [i for i, s in enumerate(result.signatures)
               if s.term[0] == sig_name]
        assert len(idx) == 1
        idx = idx[0]
    else:
        sig_name = sig_name.lower()
        idx = [i for i, s in enumerate(result.signatures)
               if s.term[3].lower().startswith(sig_name)]
        if len(idx) == 0:
            print('Error: signature name not found.', file=sys.stderr)
            return 1
        elif len(idx) > 1:
            print('Error: signature name not unique, matched: %s'
                  % ', '.join([s.term[3] for s in result.signatures]),
                  file=sys.stderr)
            return 1
        idx = idx[0]
    sig = result.signatures[idx]
    expr = result.S[idx, :]

    Y = result.Y
    y1 = Y[:, comps[0]-1]
    y2 = Y[:, comps[1]-1]

    # get signature gene expression matrix and cluster rows
    # expr = np.mean(common.get_standardized_matrix(sig.E),axis=0)

    logger.info('Plotting...')

    # plotting
    import matplotlib as mpl
    if mpl_backend is not None:
        mpl.use(mpl_backend)
    import matplotlib.pyplot as plt
    from matplotlib import rc

    if use_tex:
        rc('text', usetex=True)
    rc('font', family=font_family, size=font_size)
    rc('figure', figsize=(fig_size[0], fig_size[1]))
    rc('savefig', dpi=fig_res)

    # if fig_use_tex: rc('text', usetex = True)
    # rc('font', family = fig_font_family, size = fig_font_size)
    # rc('figure', figsize = (fig_size[0], fig_size[1]))
    # rc('savefig', dpi = fig_res)

    print(expr.shape, expr.dtype)
    print(y1.shape, y2.shape)
    # plt.scatter(y1, y2, c = expr,s=fig_marker_size,alpha=fig_marker_alpha,
    # cmap=cmap,vmin=vmin,vmax=vmax)
    plt.scatter(y1, y2, c=expr, s=200, alpha=1.0, cmap=cmap,
                vmin=vmin, vmax=vmax, edgecolor='black')

    plt.xticks(())
    plt.xlabel('PC %d Score' % comps[0])
    plt.yticks(())
    plt.ylabel('PC %d Score' % comps[1])
    plt.title(sig.get_label(include_id=False))

    # minint = int(vmin)
    # maxint = int(vmax)
    # cbticks = np.arange(int(minint),int(maxint)+0.01,1.0)
    # cb = plt.colorbar(orientation=fig_cbar_orient,shrink=fig_cbar_shrink,
    # pad=fig_cbar_pad,ticks=cbticks,use_gridspec=False,anchor=fig_cbar_anchor)
    # cb.ax.tick_params(labelsize='small')
    # cb.set_label('Standardized Expression',size='small')
    # plt.suptitle(sig_label,va='top',y=fig_title_pos)

    minint = int(vmin)
    maxint = int(vmax)
    cbticks = np.arange(minint, maxint+0.01, 1.0)
    cb = plt.colorbar(orientation=cbar_orient, shrink=cbar_scale,
                      pad=cbar_pad, ticks=cbticks, use_gridspec=False,
                      anchor=cbar_anchor)
    cb.ax.tick_params(labelsize='small')
    cb.set_label('Standardized Expression', size='small')

    simpleaxis(plt.gca())

    # plt.tight_layout()

    print('Saving to file...', end=' ')
    sys.stdout.flush()
    plt.savefig(output_file, bbox_inches='tight')
    print('done!')
    sys.stdout.flush()

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
