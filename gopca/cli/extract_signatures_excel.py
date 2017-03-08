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

"""This script extracts GO-PCA signatures as an Excel spreadsheet.

Example
-------

::

    $ gopca_extract_signatures.py -g [gopca_output_file] -o [output_file]

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
import os
# import argparse
# import cPickle as pickle
# import csv
import math

import xlsxwriter
import numpy as np

from genometools import misc
from gopca import util
from gopca.cli import arguments


def sign(x):
    return int(math.copysign(1.0, x))


def get_argument_parser():
    desc = 'Extract GO-PCA signatures as an Excel spreadsheet.'
    parser = arguments.get_argument_parser(desc=desc)
    arguments.add_io_args(parser)
    arguments.add_reporting_args(parser)
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

    # configure root logger
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose
    logger = misc.get_logger(log_file=log_file, quiet=quiet,
                             verbose=verbose)

    assert os.path.isfile(gopca_file)

    workbook = xlsxwriter.Workbook(
        output_file, {'strings_to_numbers': True, 'in_memory': True})
    workbook.set_properties({'title': 'GO-PCA Signatures'})

    bold = workbook.add_format({'bold': True})

    ws = workbook.add_worksheet()

    result = util.read_gopca_result(gopca_file)
    signatures = result.signatures

    # sort signatures first by PC, then by fold enrichment
    signatures = sorted(
        signatures, key=lambda s: [abs(s.pc), -sign(s.pc), -s.escore])

    labels = list(signatures[0].get_ordered_dict().keys())
    ws.write_row(0, 0, labels, cell_format=bold)

    max_width = np.float64([len(labels[j]) for j in range(len(labels))])
    for i, sig in enumerate(signatures):
        vals = sig.get_ordered_dict().values()
        for j, v in enumerate(vals):
            max_width[j] = max(max_width[j], float(len(v)))
        ws.write_row(i+1, 0, vals)

    for j in range(len(labels)):
        ws.set_column(j, j, max_width[j]+0.43)

    workbook.close()

    logger.info('Wrote %d signatures to "%s".', len(signatures), output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
