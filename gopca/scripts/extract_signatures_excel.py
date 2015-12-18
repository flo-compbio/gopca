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

"""This script extracts GO-PCA signatures as an Excel spreadsheet.

Example
-------

::

    $ gopca_extract_signatures.py -g [gopca_output_file] -o [output_file]

"""

import sys
import os
import argparse
import cPickle as pickle
import csv
import math

import numpy as np

from genometools import misc
from gopca import util
from gopca import cli

def get_argument_parser():
    prog = 'gopca_extract_signatures_excel.py'
    description = 'Extract GO-PCA signatures as an Excel spreadsheet.'
    parser = cli.get_argument_parser(prog, description)
    cli.add_io_args(parser)
    return parser

def main(args=None):

    sign = lambda x:int(math.copysign(1.0,x))

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    gopca_file = args.gopca_file
    output_file = args.output_file

    assert os.path.isfile(gopca_file)

    logger = misc.configure_logger(__name__)

    try:
        import xlsxwriter
    except ImportError:
        logger.error('Error: You must have the Python package "xlsxwriter" installed!')
        return 1

    workbook = xlsxwriter.Workbook(output_file,{'strings_to_numbers': True, 'in_memory': True})
    workbook.set_properties({'title': 'GO-PCA Signatures'})

    bold = workbook.add_format({'bold': True})

    ws = workbook.add_worksheet()

    output = util.read_gopca_output(gopca_file)
    signatures = output.signatures

    # sort signatures first by PC, then by fold enrichment
    signatures = sorted(signatures,
            key=lambda sig:[abs(sig.pc),-sign(sig.pc),-sig.escore])

    labels = signatures[0].get_ordered_dict().keys()
    ws.write_row(0,0,labels,cell_format=bold)

    max_width = np.float64([len(labels[j]) for j in range(len(labels))])
    for i,sig in enumerate(signatures):
        vals = sig.get_ordered_dict().values()
        for j,v in enumerate(vals):
            max_width[j] = max(max_width[j],float(len(v)))
        ws.write_row(i+1,0,vals)

    for j in range(len(labels)):
        ws.set_column(j,j,max_width[j]+0.43)

    workbook.close()

    logger.info('Wrote %d signatures to "%s".', len(signatures),output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
