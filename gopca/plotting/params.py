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

"""Functions for adding command-line parameters to GO-PCA plotting scripts.
"""

int_mv = '<int>'
float_mv = '<float>'
file_mv = '<file>'
name_mv = '<name>'

def add_fig_params(parser):
    """Add shared figure parameters.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`

    Returns
    -------
    None
    """
    g = parser.add_argument_group('Figure options')

    g.add_argument('-s', '--fig-size', type=float, nargs=2, default=[18,18],
            metavar = (float_mv, float_mv),
            help='Figure width and height (in inches).')

    g.add_argument('-r', '--fig-resolution', type=float, default=150,
            metavar = float_mv,
            help='Figure resolution (in dpi).')

    g.add_argument('-f', '--fig-font-size', type=int, default=24,
            metavar = float_mv,
            help='Figure font size (in pt).')

    g.add_argument('-m', '--fig-font-family', default='serif',
            metavar = name_mv,
            help='Figure font family.')

    g.add_argument('-t', '--fig-use-tex', action='store_true',
            help='Use LaTeX for typesetting figure text.')

    g.add_argument('-b', '--fig-mpl-backend', default=None,
            metavar = name_mv,
            help='Matplotlib backend.')

def add_heatmap_params(parser):
    """Add shared heat map parameters.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`

    Returns
    -------
    None
    """
    g = parser.add_argument_group('Heat map color and colorbar options')

    g.add_argument('-cm', '--colormap', default='RdBu_r',
            metavar = name_mv,
            help='The colormap used.')

    g.add_argument('-vc', '--val-coolest', type=float, default=-3.0,
            metavar = float_mv,
            help='The value corresponding to the "coolest" color.')

    g.add_argument('-vh', '--val-hottest', type=float, default=3.0,
            metavar = float_mv,
            help='The value corresponding to the "hottest" color.')

    g.add_argument('-co', '--cbar-orient',  default='horizontal',
            metavar = name_mv,
            help='The orientation of the colorbar.')

    g.add_argument('-ca', '--cbar-anchor', type=float, nargs=2,
            default=(0.96,1.0),
            metavar = (float_mv, float_mv),
            help='The colorbar anchoring position (x and y).')

    g.add_argument('-cs', '--cbar-scale', type=float, default=0.3,
            metavar = float_mv,
            help = 'Scaling factor to adjust the size of the colorbar.')

    g.add_argument('-cp', '--cbar-pad', type=float, default=0.015,
            metavar = float_mv,
            help = 'The colorbar padding.')


