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
import os
import argparse
import logging

import numpy as np
from sklearn.decomposition import PCA, RandomizedPCA

from genometools import misc
from gopca import common

def read_args_from_cmdline():
    parser = argparse.ArgumentParser(description='Principal Component Tester')

    parser.add_argument('-e','--expression-file',required=True)
    parser.add_argument('-G','--select-variable-genes',type=int,default=0)

    parser.add_argument('-t','--permutations',type=int,default=15)
    parser.add_argument('-z','--zscore-thresh',type=float,default=3.0)

    parser.add_argument('-s','--seed',type=int,default=None)
    parser.add_argument('--quiet',action='store_true')

    return parser.parse_args()

def main(args=None):

    if args is None:
        args = read_args_from_cmdline()

    expression_file = args.expression_file
    sel_var_genes = args.select_variable_genes
    t = args.permutations
    zscore_thresh = args.zscore_thresh
    seed = args.seed
    quiet = args.quiet

    # configure logger
    log_level = logging.INFO
    if quiet:
        log_level = logging.WARNING
    logger = misc.configure_logger(__name__, log_level = log_level)

    # set seed for random number generator
    if seed is None:
        seed = np.random.randint(int(1e9))
    np.random.seed(seed)
    logger.info('Using seed: %d', seed)

    # checks
    assert os.path.isfile(expression_file)
    assert t >= 2
    assert zscore_thresh >= 0

    # read expression
    genes,samples,E = common.read_expression(expression_file)

    # filter for most variable genes
    if sel_var_genes > 0:
        var = np.var(E,axis=1)
        a = np.argsort(var)
        a = a[::-1]
        genes = [genes[i] for i in a[:sel_var_genes]]
        E = E[a[:sel_var_genes]]
        
    logger.info('Expression matrix shape: %s', str(E.shape))
    #E += (np.random.rand(*E.shape)*1e-4)

    # do PCA on unpermuted data
    p,n = E.shape
    n_comps = min(p,n-1)
    M = PCA(n_components = n_comps)
    M.fit(E.T)
    d = M.explained_variance_ratio_.copy()

    # get permutation-based threshold
    logger.info('Performing permutations...')
    thresh = common.get_pc_explained_variance_threshold(E,zscore_thresh,t,seed)

    significant = np.sum(d >= thresh)
    args.result = significant
    logger.info('Number of significant principal components with z-score >= %.1f: %d',
            zscore_thresh, significant)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
