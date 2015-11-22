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
import cPickle as pickle

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import pearsonr

from genometools import misc
from goparser import GOParser
from gopca import common

def read_args_from_cmdline():

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-g','--gopca-file',required=True)
    parser.add_argument('-b','--bootstrap-gopca-files',required=True,nargs='+')
    parser.add_argument('-e','--expression-file',required=True)
    parser.add_argument('-t','--ontology-file',required=True)
    parser.add_argument('-o','--output-file',required=True)

    parser.add_argument('--part-of-cc-only',action='store_true')

    parser.add_argument('-fs','--figure-size',type=float,help='in inches',nargs=2,default=[18,18])
    parser.add_argument('-fr','--figure-resolution',type=int,help='in dpi',default=150)
    parser.add_argument('-ff','--figure-font-size',type=int,help='in pt',default=32)
    parser.add_argument('-fm','--figure-font-family',default='serif')
    parser.add_argument('-ft','--figure-use-tex',action='store_true')
    parser.add_argument('-fb','--figure-backend',default=None)

    parser.add_argument('-fc','--figure-colormap',default='BuPu')

    parser.add_argument('-r','--reverse-signature-order',action='store_true')
    #parser.add_argument('-fvn','--figure-vmin',type=float,default=-3.0)
    #parser.add_argument('-fvx','--figure-vmax',type=float,default=3.0)

    return parser.parse_args()

def get_signature_matrix_fast(genes,E,signatures):
    S = []
    for sig in signatures:
        indices = np.int64([misc.bisect_index(genes,g) for g in sig.genes])
        S.append(np.mean(E[indices,:],axis=0))
    S = np.float64(S)
    return S

def main(args=None):

    if args is None:
        args = read_args_from_cmdline()

    gopca_file = args.gopca_file
    bootstrap_gopca_files = args.bootstrap_gopca_files
    expression_file = args.expression_file
    ontology_file = args.ontology_file
    output_file = args.output_file
    part_of_cc_only = args.part_of_cc_only

    # figure size
    fig_size = args.figure_size
    fig_res = args.figure_resolution

    # figure text
    fig_use_tex = args.figure_use_tex
    fig_font_size = args.figure_font_size
    fig_font_family = args.figure_font_family

    # figure colormap & colorbar
    fig_cmap = args.figure_colormap

    # other
    reverse_signature_order = args.reverse_signature_order

    mpl_backend = args.figure_backend

    # read GO-PCA result
    gopca_result = common.read_gopca_result(gopca_file)
    signatures = gopca_result.signatures
    S = gopca_result.S
    #samples = gopca_result.samples

    # read expression
    genes,samples,E = common.read_expression(expression_file)
    assert np.all(np.lexsort([samples]) == np.lexsort([gopca_result.samples]))

    # sort genes + standardize expression
    a = np.lexsort([genes])
    genes_sorted = [genes[i] for i in a]
    E_sorted = common.get_standardized_matrix(E[a,:])

    # cluster signatures
    a = common.cluster_signatures(S)
    if reverse_signature_order:
        a = a[::-1]
    signatures = [signatures[i] for i in a]
    S = S[a,:]

    GO = GOParser()
    GO.parse_ontology(ontology_file,part_of_cc_only)

    # determine set of "acceptable" GO terms for each signature
    sig_term_ids = [sig.term[0] for sig in signatures]
    valid_term_ids = [set([id_]) | GO.terms[id_].ancestors | GO.terms[id_].children for id_ in sig_term_ids]
    #valid_term_ids = [set([id_]) for id_ in sig_term_ids]

    q = len(signatures)
    k = len(bootstrap_gopca_files)

    repeats = np.zeros(k,dtype=np.float64)
    resample_sizes = np.zeros(k,dtype=np.float64)

    det_med = np.zeros(k,dtype=np.float64)
    det_lq = np.zeros(k,dtype=np.float64)
    det_uq = np.zeros(k,dtype=np.float64)

    corr_med = np.zeros(k,dtype=np.float64)
    corr_lq = np.zeros(k,dtype=np.float64)
    corr_uq = np.zeros(k,dtype=np.float64)

    D = np.zeros((q,k),dtype=np.float64)
    C = np.zeros((q,k),dtype=np.float64) - 1.0

    for h,fn in enumerate(bootstrap_gopca_files):
        # for each bootstrap run

        bootstrap_result = common.read_gopca_result(fn)
        resample_sizes[h] = bootstrap_result.resample_size
        T = bootstrap_result.T
        C_max = np.zeros((T,q),dtype=np.float64) - 1.0

        found = np.zeros((T,q),dtype=np.float64)
        for j,result in enumerate(bootstrap_result.gopca_results):
            
            # generate signature matrix for signatures in bootstrapped result
            S_boot = get_signature_matrix_fast(genes_sorted,E_sorted,result.signatures)

            # calculate correlations
            S_full = np.row_stack([S,S_boot])
            C_j = np.corrcoef(S_full)[:q,q:]
            assert C_j.shape == (q,result.q)

            C_max[j,:] = np.amax(C_j,axis=1)

            for i,sig in enumerate(result.signatures):
                for i_ref in range(q):
                    if sig.term[0] in valid_term_ids[i_ref]:
                        found[j,i_ref] = 1.0

        D[:,h] = 100*np.mean(found,axis=0)
        C[:,h] = np.mean(C_max,axis=0)

    labels = [sig.get_label(include_id=False,max_name_length=50) for sig in signatures]

    import matplotlib as mpl
    if mpl_backend is not None:
        mpl.use(mpl_backend)
    import matplotlib.pyplot as plt
    from matplotlib import rc

    if fig_use_tex: rc('text',usetex=True)
    rc('font',family=fig_font_family,size=fig_font_size)
    rc('figure',figsize=(fig_size[0],fig_size[1]))
    rc('savefig',dpi=fig_res)

    plt.subplot(1,2,1)
    plt.imshow(D,interpolation='none',aspect='auto',vmin=0,vmax=100.0,cmap=fig_cmap)
    plt.yticks(np.arange(q),labels,size='x-small')
    plt.xticks(np.arange(k),['%.0f' %(s) for s in resample_sizes],size='small')
    cb = plt.colorbar(shrink=0.4)
    if fig_use_tex:
        plt.xlabel(r'\parbox{\textwidth}{\centering Sample Size\\(\% Original Sample Size)}',size='small')
        cb.set_label(r'\parbox{\textwidth}{\centering Detection Rate\\[-0.3ex] (\% Bootstrap Samples)}',size='small')
    else:
        plt.xlabel('Sample Size\n(% Original Sample Size)',size='small')
        cb.set_label('Detection Rate\n (% Bootstrap Samples)',size='small')
    #plt.title('GO Term Recovery',y=1.02)

    plt.subplot(1,2,2)
    #plt.imshow(C,interpolation='none',aspect='auto',vmin=0,vmax=1.0,cmap='Reds')
    #plt.imshow(C,interpolation='none',aspect='auto',vmin=0,vmax=1.0,cmap='jet')
    plt.imshow(C,interpolation='none',aspect='auto',vmin=0,vmax=1.0,cmap='coolwarm')
    plt.yticks(np.arange(q),())
    plt.xticks(np.arange(k),['%.0f' %(s) for s in resample_sizes],size='small')
    cb = plt.colorbar(shrink=0.4)
    if fig_use_tex:
        plt.xlabel(r'\parbox{\textwidth}{\centering Sample Size\\(\% Original Sample Size)}',size='small')
        cb.set_label(r'\parbox{\textwidth}{\centering Max. Correlation\\[-0.3ex] (Bootstrap Average)}',size='small')
    else:
        plt.xlabel('Sample Size (% Original Sample Size)',size='small')
        cb.set_label('Max. Correlation\\ (Bootstrap Average)',size='small')
    #plt.title('Signature Recovery',y=1.02)

    print 'Saving to file...', ; sys.stdout.flush()
    #plt.gcf().set_size_inches(fig_dim[0],fig_dim[1])
    plt.savefig(output_file,bbox_inches='tight')
    print 'done!'; sys.stdout.flush()
    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
