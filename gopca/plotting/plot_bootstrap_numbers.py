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

from genometools import misc
from gopca import common

def read_args_from_cmdline():

    parser = argparse.ArgumentParser(description='')

    #parser.add_argument('-g','--gopca-file',required=True)
    parser.add_argument('-b','--bootstrap-gopca-files',required=True,nargs='+')
    parser.add_argument('-o','--output-file',required=True)

    parser.add_argument('-fs','--figure-size',type=float,help='in inches',nargs=2,default=[18,12])
    parser.add_argument('-fr','--figure-resolution',type=int,help='in dpi',default=150)
    parser.add_argument('-ff','--figure-font-size',type=int,help='in pt',default=32)
    parser.add_argument('-fm','--figure-font-family',default='serif')
    parser.add_argument('-ft','--figure-use-tex',action='store_true')
    parser.add_argument('-fb','--figure-backend',default=None)
    #parser.add_argument('-fc','--figure-colormap',default='RdBu_r')
    #parser.add_argument('-fvn','--figure-vmin',type=float,default=-3.0)
    #parser.add_argument('-fvx','--figure-vmax',type=float,default=3.0)
    parser.add_argument('--pc-max',type=int,default=20)
    parser.add_argument('--sig-max',type=int,default=50)

    return parser.parse_args()

def main(args=None):

    if args is None:
        args = read_args_from_cmdline()

    #gopca_file = args.gopca_file
    bootstrap_gopca_files = args.bootstrap_gopca_files
    output_file = args.output_file

    # figure size
    fig_size = args.figure_size
    fig_res = args.figure_resolution

    # figure text
    fig_use_tex = args.figure_use_tex
    fig_font_size = args.figure_font_size
    fig_font_family = args.figure_font_family

    mpl_backend = args.figure_backend

    pc_max = args.pc_max
    sig_max = args.sig_max

    #result = common.read_gopca_result(gopca_file)

    import matplotlib as mpl
    if mpl_backend is not None:
        mpl.use(mpl_backend)
    import matplotlib.pyplot as plt
    from matplotlib import rc

    if fig_use_tex: rc('text',usetex=True)
    rc('font',family=fig_font_family,size=fig_font_size)
    rc('figure',figsize=(fig_size[0],fig_size[1]))
    rc('savefig',dpi=fig_res)

    k = len(bootstrap_gopca_files)

    resample_sizes = np.zeros(k,dtype=np.float64)
    comps_med = np.zeros(k,dtype=np.float64)
    comps_lq = np.zeros(k,dtype=np.float64)
    comps_uq = np.zeros(k,dtype=np.float64)
    sigs_med = np.zeros(k,dtype=np.float64)
    sigs_lq = np.zeros(k,dtype=np.float64)
    sigs_uq = np.zeros(k,dtype=np.float64)

    for h,fn in enumerate(bootstrap_gopca_files):
        bootstrap_result = common.read_gopca_result(fn)
        resample_sizes[h] = bootstrap_result.resample_size
        T = bootstrap_result.T
        comps = np.zeros(T,dtype=np.float64)
        sigs = np.zeros(T,dtype=np.float64)
        for j,result in enumerate(bootstrap_result.gopca_results):
            comps[j] = result.D
            sigs[j] = result.q
        comps_med[h] = np.median(comps)
        comps_lq[h] = np.percentile(comps,25.0)
        comps_uq[h] = np.percentile(comps,75.0)
        sigs_med[h] = np.median(sigs)
        sigs_lq[h] = np.percentile(sigs,25.0)
        sigs_uq[h] = np.percentile(sigs,75.0)

    print comps_lq
    print comps_med
    print comps_uq
    error_kw = {'ecolor': 'gray', 'lw': 3.0, 'capsize': 10, 'mew': 3.0, 'zorder': 50}
    plt.bar(np.arange(k)-0.4,comps_med,yerr=[comps_med-comps_lq,comps_uq-comps_med],\
            width=0.8,edgecolor='none',color='skyblue',error_kw=error_kw)
    plt.ylim(0,pc_max)
    if fig_use_tex:
        plt.ylabel(r'\# Principal Components Tested')
    else:
        plt.ylabel('# Principal Components Tested')
    common.simpleaxis(plt.gca())

    if fig_use_tex:
        plt.xlabel(r'Sample Size (\% Original Sample Size)')
    else:
        plt.xlabel('Sample Size (% Original Sample Size)')

    ax2 = plt.twinx()
    plt.sca(ax2)

    error_kw = {'ecolor': 'purple', 'lw': 3.0, 'capsize': 10, 'mew': 3.0, 'zorder': 50}
    plt.errorbar(np.arange(k),sigs_med,yerr=[sigs_med-sigs_lq,sigs_uq-sigs_med],\
            color='purple',**error_kw)

    plt.xticks(np.arange(k),['%.0f' %(s) for s in resample_sizes])
    plt.xlim(-0.5,k-0.5)
    if fig_use_tex:
        plt.ylabel(r'\# Signatures Generated',color='purple')
    else:
        plt.ylabel('# Signatures Generated',color='purple')
    plt.yticks(color='purple')
    plt.ylim(0,sig_max)
    plt.gca().spines['top'].set_visible(False)

    #for r in range(T):
    #   comps[r] = np.median(C[r,:])
    #   lq[h] = comps[h] - np.percentile(C[h,:],25.0)
    #   uq[h] = np.percentile(C[h,:],75.0) - comps[h]
    print 'Saving to file...', ; sys.stdout.flush()
    #plt.gcf().set_size_inches(fig_dim[0],fig_dim[1])
    plt.savefig(output_file,bbox_inches='tight')
    print 'done!'; sys.stdout.flush()
    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
