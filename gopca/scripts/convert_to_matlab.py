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

import numpy as np
from scipy.io import savemat

from gopca import util
from gopca import cli

def get_argument_parser():

    prog = 'convert_to_matlab_format.py'
    description = 'Converts GO-PCA output to MATLAB format.'
    parser = cli.get_argument_parser(prog, description)

    cli.add_io_args(parser)
    #params.add_reporting_params(parser)

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
    #log_file = args.log_file
    #quiet = args.quiet
    #verbose = args.verbose

    G = util.read_gopca_output(gopca_file)
    
    signatures = G.signatures
    #for sig in signatures:
    #    sig.genes = np.asarray(sig.genes, dtype = np.object)
    #    sig.enr.genes = np.asarray(sig.enr.genes, dtype = np.object)
    #    sig.enr.term = np.asarray(sig.enr.term, dtype = np.object)

    sig_labels = np.asarray([sig.label for sig in signatures],
            dtype = np.object)
    sig_genes = [np.asarray(sig.genes, dtype = np.object) for sig in signatures]
    sig_term_genes = [np.asarray(sig.enr.genes, dtype = np.object) for sig in signatures]
    #sig_genes = dict([sig.term[0].replace(':','_'),sorted(sig.genes)] for sig in signatures)
    #print sig_genes

    #common.write_expression(output_file,labels,samples,S)

    mat = {}

    # configuration (savemat does not know how to handle None)
    conf = G.user_config.get_dict()
    for k,v in conf.iteritems():
        if v is None:
            conf[k] = ''
    mat['user_config'] = conf

    conf = G.config.get_dict()
    for k,v in conf.iteritems():
        if v is None:
            conf[k] = ''
    mat['config'] = conf

    # gene and sample names
    mat['genes'] = np.asarray(G.genes, dtype = np.object)
    mat['samples'] = np.asarray(G.samples, dtype = np.object)

    # PCA data
    mat['W'] = G.W # the PC loading matrix
    mat['Y'] = G.Y # the PC score matrix

    # signatures and signature matrix
    mat['signatures'] = signatures
    mat['signature_labels'] = sig_labels
    mat['signature_genes'] = sig_genes
    mat['signature_term_genes'] = sig_term_genes

    mat['signature_genes'] = sig_genes
    mat['S'] = G.S # the signature matrix

    # output hash
    mat['output_hash'] = G.hash

    # signatures
    savemat(output_file, mat, appendmat = append_mat)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
