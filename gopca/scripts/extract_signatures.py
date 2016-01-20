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

"""This script extracts GO-PCA signatures as a tab-delimited text file.

Example
-------

::

    $ gopca_extract_signatures.py -g [gopca_output_file] -o [output_file]

"""

import sys
import os
import argparse
import csv
import math

import numpy as np

from genometools import misc
from gopca import util
from gopca import cli

sign = lambda x:int(math.copysign(1.0, x))

def get_argument_parser():

    desc = 'Extract GO-PCA signatures as a tab-delimited text file.'
    parser = cli.get_argument_parser(desc = desc)
    cli.add_io_args(parser)

    return parser

def main(args=None):

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    gopca_file = args.gopca_file
    output_file = args.output_file

    logger = misc.configure_logger(__name__)

    assert os.path.isfile(gopca_file)

    G = util.read_gopca_result(gopca_file)
    signatures = G.signatures

    # sort signatures first by PC, then by fold enrichment
    signatures = sorted(signatures, 
            key=lambda sig:[abs(sig.pc),-sign(sig.pc),-sig.escore])

    labels = signatures[0].get_ordered_dict().keys()

    with open(output_file,'w') as ofh:
        writer = csv.writer(ofh,dialect='excel-tab',lineterminator='\n',
                quoting=csv.QUOTE_NONE)

        writer.writerow(labels)

        for i,sig in enumerate(signatures):
            vals = sig.get_ordered_dict().values()
            writer.writerow(vals)

    logger.info('Wrote %d signatures to "%s".',
            len(signatures), output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
