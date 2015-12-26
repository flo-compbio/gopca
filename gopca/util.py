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

"""Functions used by various GO-PCA scripts.
"""

import os
import sys
import argparse
import cPickle as pickle
import hashlib
import logging
from pkg_resources import parse_version

import unicodecsv as csv

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import sklearn

import gopca
from genometools import misc

#logger = logging.getLogger(__name__)

def get_logger(name = '', log_file = None, quiet = False,
    verbose = False):

    # configure root logger
    log_level = logging.INFO
    if quiet:
        log_level = logging.WARNING
    elif verbose:
        log_level = logging.DEBUG

    new_logger = misc.configure_logger(name, log_file = log_file,
            log_level = log_level)

    return new_logger

def filter_signatures(signatures, S, corr_thresh):

    if corr_thresh == 1.0:
        # no filtering
        return signatures, S

    q, n = S.shape
    info = np.zeros(q, dtype = np.float64)
    # sort signatures by information content?
    max_ent = - n * ((1/float(n)) * np.log2(1/float(n)))
    for i, sig in enumerate(signatures):
        abs_sig = np.absolute(S[i,:])
        total = np.sum(abs_sig)
        rel = abs_sig / total
        ent = - np.sum(rel * np.log2(rel))
        info[i] = max_ent - ent
    a = np.argsort(info)
    a = a[::-1]

    sig_abs_pcs = np.absolute(np.int64([sig.pc for sig in signatures]))
    a = np.lexsort([-info, sig_abs_pcs])
    #print '\n'.join(['%.3f: %s' %(info[i],str(G.signatures[i])) for i in a[:5]])

    # filtering
    sel = np.ones(q, dtype = np.bool_)
    for i in a:

        if not sel[i]:
            # already excluded
            continue

        for i2, sig in enumerate(signatures):
            if i == i2 or not sel[i2]:
                continue
            assert np.corrcoef(np.vstack([S[i,:],S[i2,:]])).shape == (2,2)
            if np.corrcoef(np.vstack([S[i,:],S[i2,:]]))[0,1] >= corr_thresh:
                sel[i2] = False

    sel = np.nonzero(sel)[0]
    signatures = [signatures[i] for i in sel]
    S = S[sel,:]

    return signatures, S

def get_pc_explained_variance_threshold(E, z, t, seed):

    # RandomizedPCA does not work in Scikit-learn 0.14.1,
    # but it works in Scikit-learn 0.16.1
    if parse_version(sklearn.__version__) >= parse_version('0.16.1'):
        from sklearn.decomposition import RandomizedPCA as PCA
    else:
        from sklearn.decomposition import PCA

    # initialize random number generator
    np.random.seed(seed)

    # do permutations
    p,n = E.shape
    d_max_null = np.empty(t,dtype=np.float64)
    E_perm = np.empty((p,n),dtype=np.float64)
    M_null = PCA(n_components = 1)
    for j in xrange(t):

        for i in xrange(p):
            E_perm[i,:] = E[i,np.random.permutation(n)]

        M_null.fit(E_perm.T)
        d_max_null[j] = M_null.explained_variance_ratio_[0]

    # calculate z-score threshold
    mean_null = np.mean(d_max_null)
    std_null = np.std(d_max_null,ddof=1)
    thresh = mean_null + z * std_null

    return thresh

def get_file_md5sum(path, mode = 'rb'):
    """Get MD5 hash of file content.

    Parameters
    ----------
    path: str
        Path of file.
    
    Returns
    -------
    str
        MD5 hash of file content, represented as a 32-digit hex string.
    """
    digest = None
    with open(path, mode = mode) as fh:
        digest = hashlib.md5(fh.read()).hexdigest()
    return digest
    
def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def print_signatures(signatures):
    a = None
    maxlength = 40
    a = sorted(range(len(signatures)),key=lambda i: -signatures[i].msfe)

    for i in a:
        sig = signatures[i]
        print sig.get_label(max_name_length=maxlength,include_pval=True)

def get_centered(e):
    e = e.copy()
    e -= np.mean(e)
    return e

def get_centered_matrix(E):
    return np.float64([get_centered(e) for e in E])

def get_standardized(e):
    e = e.copy()
    e -= np.mean(e)
    e /= np.std(e,ddof=1)
    return e

def get_standardized_matrix(E):
    return np.float64([get_standardized(e) for e in E])

#def get_mean_standardized_(E):
#   return np.mean(np.float64([get_standardized(e) for e in E]),axis=0)

def get_signature_expression(genes,E,sig_genes):
    p_sig = len(sig_genes)
    p,n = E.shape
    S = np.zeros((p_sig,n),dtype=np.float64)
    for i,g in enumerate(sig_genes):
        idx = genes.index(g)
        S[i,:] = E[idx,:]
        S[i,:] -= np.mean(S[i,:])
        S[i,:] /= np.std(S[i,:],ddof=1)
    sig = np.mean(S,axis=0)
    return sig

def get_signature_expression_robust(genes,E,sig_genes):
    p_sig = len(sig_genes)
    p,n = E.shape
    S = np.zeros((p_sig,n),dtype=np.float64)
    for i,g in enumerate(sig_genes):
        idx = misc.bisect_index(genes,g)
        S[i,:] = E[idx,:]
        med = np.median(S[i,:])
        mad = np.median(np.absolute(S[i,:]-med))
        std = 1.4826*mad
        S[i,:] -= med
        S[i,:] /= std
    sig = np.mean(S,axis=0)
    return sig

def get_median_pairwise_correlation(E):
    C = np.corrcoef(E)
    sel = np.triu_indices(C.shape[0], k=1)
    return np.median(C[sel])

def get_signature_label(GO, sig, max_length=40):
    count = ' (%d:%d/%d)' %(sig.pc,len(sig.genes),sig.K)
    enr = sig.enr
    return GO.terms[enr.term[0]].get_pretty_format(omit_acc=True,max_name_length=max_length) + count

def variance_filter(genes, E, top):
    # filter genes by variance
    a = np.argsort(np.var(E,axis=1,ddof=1))[::-1]
    n = E.shape[0]
    sel = np.zeros(n,dtype=np.bool_)
    sel[a[:top]] = True
    sel = np.nonzero(sel)[0]
    genes = [genes[i] for i in sel]
    E = E[sel,:]
    return genes,E

def read_gopca_output(path):
    """Read GO-PCA output from pickle."""
    G = None
    with open(path, 'rb') as fh:
        G = pickle.load(fh)
    if isinstance(G, gopca.GOPCARun):
        G = G.output
    assert isinstance(G, gopca.GOPCAOutput)
    return G

def read_go_annotations(fn):
    ann = {}
    with open(fn, 'rb') as fh:
        reader = csv.reader(fh, dialect='excel-tab')
        for l in reader:
            ann[tuple(l[:4])] = l[4].split(',')
    return ann

def cluster_rows(S, metric='correlation', method='average', reverse=False):
    distxy = squareform(pdist(S, metric=metric))
    R = dendrogram(linkage(distxy, method=method), no_plot=True)
    order_rows = np.int64([int(l) for l in R['ivl']])
    if reverse:
        order_rows = order_rows[::-1]
    return order_rows

def cluster_signatures(S, metric = 'correlation', method = 'average',
        reverse = False):
    # hierarchical clustering of signatures
    order_rows = cluster_rows(S, metric, method, reverse)
    return order_rows

#def cluster_samples(S, metric = 'euclidean', method = 'average',
def cluster_samples(S, metric = 'correlation', method = 'average',
        reverse = False):
    order_cols = cluster_rows(S.T, metric, method, reverse)
    return order_cols

def get_qvalues(pvals, pi_zero = 1.0):
    # implements storey-tibshirani procedure for calculating q-values
    n = pvals.size
    qvals = np.empty(n,dtype=np.float64)

    # start with largest p-value
    a = np.argsort(pvals,kind='mergesort') # stable sort
    a = a[::-1]

    s = 1
    q = 1.0
    for i in a:
        q = min(((pi_zero * pvals[i])*n)/s , q)
        qvals[i] = q
        s += 1

    return qvals
