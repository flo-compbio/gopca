#!/usr/bin/env python

# Copyright (c) 2017 Florian Wagner
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

"""This script generates an interactive plot of the GO-PCA signature matrix.

The output format is HTML! A png file can be downloaded manually by opening
the HTML file in a browser and clicking on the cameral symbol in the top right
corner of the figure.

Example
-------

::

    $ gopca_plot_signature_matrix.py -g gopca_result.pickle \
            -o gopca_sig_matrix.html

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import sys

from plotly.offline import plot

from genometools import misc
from gopca import util
from gopca.cli import arguments




def get_argument_parser():
    desc = 'Generate an interactive plot of the GO-PCA signature matrix.'
    parser = arguments.get_argument_parser(desc=desc)
    arguments.add_io_args(parser)
    arguments.add_reporting_args(parser)
    # arguments.add_sample_args(parser)
    arguments.add_figure_args(parser)
    arguments.add_heatmap_args(parser)

    g = parser.add_argument_group('Plotting options')
    g.add_argument(
        '--no-plotly-js', action='store_true',
        help='Do not include plotly javascript code in figure.')
    
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

    emin = args.min_val
    emax = args.max_val

    width = args.width
    height = args.height

    margin_left = args.margin_left
    margin_bottom = args.margin_bottom

    font_size = args.font_size
    font = args.font

    include_plotlyjs = True
    if args.no_plotly_js:
        include_plotlyjs = False

    # configure root logger
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose
    logger = misc.get_logger(log_file=log_file, quiet=quiet,
                             verbose=verbose)

    colorbar_label = args.colorbar_label

    sig_matrix = util.read_gopca_result(gopca_file)
    fig = sig_matrix.get_figure(
        width=width, height=height,
        font_size=font_size, font=font,
        emin=emin, emax=emax,
        margin_left=margin_left, margin_bottom=margin_bottom,
        heatmap_kw=dict(colorbar_label=colorbar_label))
    plot(fig, filename=output_file, image_filename='gopca_signature_matrix', 
         auto_open=False, include_plotlyjs=include_plotlyjs)

    logger.info('Plotted  %d x %d signature matrix.',
                sig_matrix.p, sig_matrix.n)

    return 0


if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)