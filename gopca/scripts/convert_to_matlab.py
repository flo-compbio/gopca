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

import sys
import argparse
import cPickle as pickle

from scipy.io import savemat

from gopca import util
from gopca import params

def get_argument_parser():

    prog = 'convert_to_matlab_format.py'
    description = 'Converts GO-PCA output to MATLAB format.'
    parser = params.get_argument_parser(prog, description)

    params.add_io_params(parser)
    params.add_reporting_params(parser)

    parser.add_argument('--append-mat', action = 'store_true',
            help = 'Automatically append .mat file extension.')

    return parser

def main(args=None):

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    # input and output files
    gopca_file = args.gopca_file
    output_file = args.output_file

    # specific options
    append_mat = args.append_mat

    # reporting options
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    G = util.read_gocpca_output(gopca_file)
    
    signatures = G.signatures
    labels = [sig.get_label() for sig in signatures]
    sig_genes = dict([sig.term[0].replace(':','_'),sorted(sig.genes)] for sig in signatures)
    #print sig_genes
    samples = list(result.samples)
    S = result.S
    #common.write_expression(output_file,labels,samples,S)

    mat = {}

    # parameters, file hashes and global hash
    mat['params'] = G.get_input_dict()
    mat['raw_params'] = G.get_raw_input_dict()

    # gene and sample names
    mat['genes'] = G.genes
    mat['samples'] = G.samples

    # PCA data
    mat['W'] = G.W # the PC loading matrix
    mat['Y'] = G.Y # the PC score matrix

    # signatures and signature matrix
    mat['signature_labels'] = labels
    mat['signature_genes'] = sig_genes

    mat['S'] = G.S # the signature matrix

    # signatures
    savemat(output_file,mat,appendmat=append_mat)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
