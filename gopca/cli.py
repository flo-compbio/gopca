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

import sys
import argparse
import textwrap

import gopca

file_mv = '<file>'
int_mv = '<int>'
float_mv = '<float>'
name_mv = '<name>'
str_mv = '<str>'

str_type = lambda s: unicode(s, sys.getfilesystemencoding())

def get_argument_parser(prog = None, desc = None,
        formatter_class = None):
    """Create an argument parser.

    Parameters
    ----------
    prog: str, optional
        The program name.
    description: str, optinoal
        The program description.

    Returns
    -------
    `argparse.ArgumentParser`
        The arguemnt parser created.
    """
    if formatter_class is None:
        formatter_class = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(prog = prog, description = desc,
            formatter_class = formatter_class, add_help = False)

    g = parser.add_argument_group('Help')
    g.add_argument('-h', '--help', action='help',
            help = 'Show this help message and exit.')

    v = gopca.__version__
    g.add_argument('--version', action='version', version='GO-PCA ' + v,
            help = 'Output the GO-PCA version and exit.')

    return parser

def add_reporting_args(parser):
    """Add reporting arguments to an argument parser.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`

    Returns
    -------
    `argparse.ArgumentGroup`
        The argument group created.
    """
    g = parser.add_argument_group('Reporting options')

    g.add_argument('-l', '--log-file', default = None,
            metavar = file_mv, type = str_type, help = textwrap.dedent("""\
                Path of log file (if specified, report to stdout AND
                file)."""))

    g.add_argument('-q', '--quiet', action = 'store_true',
            help = 'Only output errors and warnings.')

    g.add_argument('-v', '--verbose', action = 'store_true',
            help = 'Enable verbose output. Ignored if --quiet is specified.')

    return parser

def add_io_args(parser):
    """Add input/output arguments to an argument parser.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`

    Returns
    -------
    `argparse.ArgumentGroup`
        The argument group created.
    """
    g = parser.add_argument_group('Input and output files (required)')

    g.add_argument('-g', '--gopca-file', required = True,
            metavar = file_mv, type = str_type,
            help = 'The GO-PCA output file.')

    g.add_argument('-o', '--output-file', required = True,
            metavar = file_mv, type = str_type,
            help = 'The output file.')

    return g

def add_go_term_args(parser):
    """Add GO term arguments to an argument parser.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`

    Returns
    -------
    `argparse.ArgumentGroup`
        The argument group created.
    """
    g = parser.add_argument_group('GO term options')

    g.add_argument('--term-reverse-order', action = 'store_true',
            help='Reverse the order of the GO terms.')

    g.add_argument('--term-max-len', type = int, default = 50,
            metavar = int_mv,
            help='The maximal length of GO term labels.')

    return g

def add_signature_args(parser):
    """Add signature arguments to an argument parser.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`

    Returns
    -------
    `argparse.ArgumentGroup`
        The argument group created.
    """
    g = parser.add_argument_group('Signature options')

    g.add_argument('--sig-reverse-order', action = 'store_true',
            help='Reverse the order of the signatures.')

    g.add_argument('--sig-max-len', type = int, default = 50,
            metavar = int_mv,
            help='The maximal length of signature labels.')

    g.add_argument('--sig-filter-corr', type = float, default = 1.0,
            metavar = float_mv,
            help = textwrap.dedent("""\
                Correlation threshold for filtering signatures
                (1.0 = off)."""))

    return g

def add_sample_args(parser):
    """Add sample arguments to an argument parser.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`

    Returns
    -------
    `argparse.ArgumentGroup`
        The argument group created.
    """
    g = parser.add_argument_group('Sample options')

    g.add_argument('--sample-no-clustering', action = 'store_true',
            help = 'Disable clustering of the samples.')

    g.add_argument('--sample-cluster-metric', default='correlation',
            metavar = name_mv, help = textwrap.dedent("""\
                The metric used in the hierarchical clustering algorithm."""))

    return g
