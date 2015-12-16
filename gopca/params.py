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

def add_reporting_params(parser):
    # reporting options
    g = parser.add_argument_group('Reporting options')

    g.add_argument('-l', '--log-file', default=None,
            metavar = file_mv,
            help = 'Path of log file (if specified, report to stdout AND ' +
            'file.')

    g.add_argument('-q', '--quiet', action='store_true',
            help = 'Only output errors and warnings.')

    g.add_argument('-v', '--verbose', action='store_true',
            help = 'Enable verbose output. Ignored if --quiet is specified.')

    return parser

def add_io_params(parser):
    """Add shared input/output parameters to an argument parser.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`

    Returns
    -------
    None
    """
    g = parser.add_argument_group('Input and output files (required)')

    g.add_argument('-g', '--gopca-file', required=True,
            metavar = file_mv,
            help = 'The GO-PCA output file.')

    g.add_argument('-o', '--output-file', required=True,
            metavar = file_mv,
            help = 'The output file.')

def add_go_term_params(parser):
    """Add shared go term parameters.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`

    Returns
    -------
    None
    """
    g = parser.add_argument_group('GO term options')

    g.add_argument('--term-reverse-order', action='store_true',
            help='Reverse the order of the GO terms.')

    g.add_argument('--term-max-len', type=int, default=50,
            metavar = int_mv,
            help='The maximal length of GO term labels.')

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

    g.add_argument('--sig-max-len', type=int, default=50,
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

