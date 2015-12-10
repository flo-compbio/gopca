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

"""Functions for configuring command-line parameters of GO-PCA scripts.
"""

import argparse

import gopca

file_mv = '<file>'
int_mv = '<int>'
float_mv = '<float>'
name_mv = '<name>'

def get_argument_parser(prog, description, formatter_class = None):

    if formatter_class is None:
        formatter_class = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(prog = prog, description = description,
            formatter_class = formatter_class, add_help = False)

    g = parser.add_argument_group('Help')
    g.add_argument('-h', '--help', action='help',
            help='Show this help message and exit.')

    v = gopca.__version__
    g.add_argument('--version', action='version', version='GO-PCA ' + v,
            help='Output the GO-PCA version and exit.')

    return parser

def add_signature_params(parser):
    """Add shared signature parameters.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`

    Returns
    -------
    None
    """
    g = parser.add_argument_group('Signature options')

    g.add_argument('--sig-reverse-order', action='store_true',
            help='Reverse the order of the signatures.')

    g.add_argument('-l', '--sig-max-len', type=int, default=50,
            metavar = int_mv,
            help='The maximal length of signature labels.')


def add_sample_params(parser):
    """Add shared sample parameters.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`

    Returns
    -------
    None
    """
    g = parser.add_argument_group('Sample options')

    g.add_argument('--sample-no-clustering', action='store_true',
            help='Disable clustering of the samples.')

    g.add_argument('--sample-cluster-metric', default='euclidean',
            metavar = name_mv,
            help='The metric used in the hierarchical clustering algorithm.')

