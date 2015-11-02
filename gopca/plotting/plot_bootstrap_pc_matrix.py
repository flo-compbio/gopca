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
    parser.add_argument('-b','--bootstrap-gopca-file',required=True)
    parser.add_argument('-e','--expression-file',required=True)
    parser.add_argument('-t','--ontology-file',required=True)
    parser.add_argument('-o','--output-file',required=True)

    #parser.add_argument('--identical-terms-only',action='store_true') # to-do!

    parser.add_argument('-D','--principal-components',type=int,required=True)
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
    bootstrap_gopca_file = args.bootstrap_gopca_file
    expression_file = args.expression_file
    ontology_file = args.ontology_file
    output_file = args.output_file
    part_of_cc_only = args.part_of_cc_only
    n_comps = args.principal_components

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
    #samples = gopca_result.samples
    signatures = gopca_result.signatures
    S = gopca_result.S

    # read bootstrapping result
    bootstrap_result = common.read_gopca_result(bootstrap_gopca_file)
    I = bootstrap_result.I

    # read expression
    genes,samples,E = common.read_expression(expression_file)
    assert np.all(np.lexsort([samples]) == np.lexsort([gopca_result.samples]))

    # sorting
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

    #repeats = np.zeros(k,dtype=np.float64)
    q = len(signatures)
    D = np.empty((q,n_comps),dtype=np.float64)
    C_max = np.empty((q,n_comps),dtype=np.float64)

    # calculate correlations first
    print 'Calculating correlations...', ; sys.stdout.flush()
    T = bootstrap_result.T
    C = np.zeros((T,n_comps,q),dtype=np.float64) - 1.0
    for j,result in enumerate(bootstrap_result.gopca_results):
        assert result.n == gopca_result.n

        # generate signature matrix for signatures in bootstrapped result
        S_boot = get_signature_matrix_fast(genes_sorted,E_sorted,result.signatures)

        # calculate correlations
        S_full = np.row_stack([S,S_boot])
        C_j = np.corrcoef(S_full)[:q,q:]
        assert C_j.shape == (q,result.q)

        # determine maximum correlations by PC
        for i,sig in enumerate(result.signatures):
            pci = abs(sig.pc)-1
            if pci < n_comps: # for some bootstrap samples, signatures are generated for PCs larger than the number of PCs tested in the original analysis
                C[j,pci,:] = np.amax(np.c_[C[j,pci,:],C_j[:,i]],axis=1)
    print 'done!'; sys.stdout.flush()
    

    for max_pc_index in range(n_comps):
        T = bootstrap_result.T
        found = np.zeros((T,q),dtype=np.float64)
        maxcorr = np.zeros((T,q),dtype=np.float64)-1.0
        for j,result in enumerate(bootstrap_result.gopca_results):
            # for each bootstrap sample
            S_sub = S[:,I[j,:]]
            found_repeat = np.zeros(q,dtype=np.float64)
            
            for i,sig in enumerate(result.signatures):
                if abs(sig.pc)-1 <= max_pc_index:
                    for i_ref in range(q):
                        if sig.term[0] in valid_term_ids[i_ref]:
                            found_repeat[i_ref] = 1.0
            found[j,:] = found_repeat
        D[:,max_pc_index] = 100*np.mean(found,axis=0) # convert to %
        maxcorr = np.amax(C[:,:(max_pc_index+1),:],axis=1)
        assert maxcorr.shape == (T,q)
        C_max[:,max_pc_index] = np.mean(maxcorr,axis=0)

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
    extent = (0.5,n_comps+0.5,q-0.5,-0.5)
    plt.imshow(D,extent=extent,interpolation='none',aspect='auto',vmin=0,vmax=100.0,cmap=fig_cmap)
    #plt.imshow(D,interpolation='none',aspect='auto',vmin=0,vmax=100.0,cmap=fig_cmap)
    xticks,_ = plt.xticks()
    xticklabels = [str(int(x)) for x in xticks]
    plt.xticks(xticks,xticklabels,size='small')
    plt.xlim(0.5,n_comps+0.5)
    plt.yticks(np.arange(q),labels,size='x-small')
    cb = plt.colorbar(shrink=0.4)
    if fig_use_tex:
        plt.xlabel(r'\parbox{\textwidth}{\centering \# Principal Components\\ Included}',size='small')
        cb.set_label(r'\parbox{\textwidth}{\centering Detection Rate\\[-0.3ex] (\% Bootstrap Samples)}',size='small')
    else:
        plt.xlabel('# Principal Components\nIncluded',size='small')
        cb.set_label('Detection Rate\n (% Bootstrap Samples)',size='small')
    cb.ax.tick_params(labelsize='small')

    plt.subplot(1,2,2)
    extent = (0.5,n_comps+0.5,q-0.5,-0.5)
    plt.imshow(C_max,extent=extent,interpolation='none',aspect='auto',vmin=0,vmax=1.0,cmap='coolwarm')
    #plt.imshow(C_max,interpolation='none',aspect='auto',vmin=0,vmax=1.0,cmap='coolwarm')
    xticks,_ = plt.xticks()
    xticklabels = [str(int(x)) for x in xticks]
    plt.xticks(xticks,xticklabels,size='small')
    plt.xlim(0.5,n_comps+0.5)
    plt.yticks(np.arange(q),())
    cb = plt.colorbar(shrink=0.4)
    if fig_use_tex:
        plt.xlabel(r'\parbox{\textwidth}{\centering \# Principal Components\\ Included}',size='small')
        cb.set_label(r'\parbox{\textwidth}{\centering Max. Correlation\\[-0.3ex] (Bootstrap Average)}',size='small')
    else:
        plt.xlabel('# Principal Components\nIncluded',size='small')
        cb.set_label('Max. Correlation\\ (Bootstrap Average)',size='small')
    cb.ax.tick_params(labelsize='small')

    print 'Saving to file...', ; sys.stdout.flush()
    #plt.gcf().set_size_inches(fig_dim[0],fig_dim[1])
    plt.savefig(output_file,bbox_inches='tight')
    print 'done!'; sys.stdout.flush()
    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
