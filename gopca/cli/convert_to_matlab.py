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

"""Script to convert GO-PCA result to matlab file format."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys

import numpy as np
from scipy.io import savemat

from gopca import util
from gopca.cli import arguments


def get_argument_parser():

    desc = 'Converts GO-PCA result to MATLAB format.'
    parser = arguments.get_argument_parser(desc=desc)

    arguments.add_io_args(parser)
    # arguments.add_reporting_args(parser)

    parser.add_argument('--append-mat', action='store_true',
                        help='Automatically append .mat file extension.')

    return parser


def main(args=None):

    vinfo = sys.version_info
    if not (vinfo >= (2, 7)):
        raise SystemError('Python interpreter version >= 2.7 required, '
                          'found %d.%d instead.' % (vinfo.major, vinfo.minor))

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    # input and output files
    gopca_file = args.gopca_file
    output_file = args.output_file

    # specific options
    append_mat = args.append_mat

    # reporting options
    # log_file = args.log_file
    # quiet = args.quiet
    # verbose = args.verbose

    res = util.read_gopca_result(gopca_file)
    
    signatures = res.signatures
    # for sig in signatures:
    #    sig.genes = np.asarray(sig.genes, dtype = np.object)
    #    sig.enr.genes = np.asarray(sig.enr.genes, dtype = np.object)
    #    sig.enr.term = np.asarray(sig.enr.term, dtype = np.object)

    sig_labels = np.asarray([sig.label for sig in signatures],
                            dtype=np.object)
    sig_genes = [np.asarray(sig.genes, dtype=np.object)
                 for sig in signatures]
    sig_term_genes = [np.asarray(sig.enr.genes, dtype=np.object)
                      for sig in signatures]
    # sig_genes = dict([sig.term[0].replace(':','_'),sorted(sig.genes)]
    #                   for sig in signatures)
    # print(sig_genes)

    # common.write_expression(output_file,labels,samples,S)

    mat = {}

    # configuration (savemat does not know how to handle None)
    conf = res.config.get_dict()
    for k, v in conf.items():
        if v is None:
            conf[k] = ''
    mat['config'] = conf

    # fix signature gene sets (savemat does not handle frozenset)
    for s in res.signatures:
        s.enr.gene_set.genes = tuple(sorted(s.enr.gene_set.genes))

    # gene and sample names
    mat['genes'] = np.asarray(res.genes, dtype=np.object)
    mat['samples'] = np.asarray(res.samples, dtype=np.object)

    # PCA data
    mat['W'] = res.W  # the PC loading matrix
    mat['Y'] = res.Y  # the PC score matrix

    # signatures and signature matrix
    mat['signatures'] = signatures
    mat['signature_labels'] = sig_labels
    mat['signature_genes'] = sig_genes
    mat['signature_term_genes'] = sig_term_genes

    mat['signature_genes'] = sig_genes
    mat['S'] = res.S  # the signature matrix

    # output hash
    mat['output_hash'] = res.hash

    # signatures
    savemat(output_file, mat, appendmat=append_mat)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
