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

"""Functions used by various GO-PCA scripts.
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

# import os
import io
import sys
# import argparse
import copy
import hashlib
import logging
# from collections import Iterable

import six
import unicodecsv as csv

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import sklearn

from genometools import misc
from genometools.expression import ExpMatrix, cluster
import gopca
from gopca import GOPCASignatureMatrix
from gopca import GOPCASignature

if six.PY2:
    import cPickle as pickle
else:
    import pickle

logger = logging.getLogger(__name__)


def combine_signatures(*args):
    """Combines signatures from multiple GO-PCA signature matrices.

    """
    # TODO: Finish docstring
    for result in args:
        assert isinstance(result, GOPCASignatureMatrix)
    # assert len(results
    merged = copy.deepcopy(args[0])
    for i, other in enumerate(args[1:]):
        if other.samples != merged.samples:
            raise ValueError('Cannot combine results with different '
                             'samples.')
        merged.signatures.extend(other.signatures)
    return merged


def get_logger(name='', log_stream=sys.stdout, log_file=None,
               quiet=False, verbose=False):

    # configure root logger
    log_level = logging.INFO
    if quiet:
        log_level = logging.WARNING
    elif verbose:
        log_level = logging.DEBUG

    new_logger = misc.configure_logger(name, log_stream, log_file, log_level)

    return new_logger




def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


def print_signatures(signatures):
    # a = None
    maxlength = 40
    a = sorted(range(len(signatures)), key=lambda i: -signatures[i].escore)

    for i in a:
        sig = signatures[i]
        print(sig.get_label(max_name_length=maxlength, include_pval=True))


def get_centered(x):
    x = x.copy()
    x -= np.mean(x)
    return x


def get_centered_matrix(X):
    return np.float64([get_centered(x) for x in X])


def get_standardized(x):
    x = x.copy()
    x -= np.mean(x)
    x /= np.std(x, ddof=1)
    return x


def get_standardized_matrix(X):
    return np.float64([get_standardized(x) for x in X])

# def get_mean_standardized_(E):
#   return np.mean(np.float64([get_standardized(e) for e in E]),axis=0)


def get_signature_expression(genes, X, sig_genes, standardize=True):
    p_sig = len(sig_genes)
    p, n = X.shape
    S = np.zeros((p_sig, n), dtype=np.float64)
    for i, g in enumerate(sig_genes):
        idx = genes.index(g)
        S[i, :] = X[idx, :]
        S[i, :] -= np.mean(S[i, :])
        if standardize:
            S[i, :] /= np.std(S[i, :], ddof=1)
    sig = np.mean(S, axis=0)
    return sig


def get_signature_expression_robust(genes, E, sig_genes):
    p_sig = len(sig_genes)
    p, n = E.shape
    S = np.zeros((p_sig, n), dtype=np.float64)
    for i, g in enumerate(sig_genes):
        idx = misc.bisect_index(genes, g)
        S[i, :] = E[idx, :]
        med = np.median(S[i, :])
        mad = np.median(np.absolute(S[i, :] - med))
        std = 1.4826*mad
        S[i, :] -= med
        S[i, :] /= std
    sig = np.mean(S, axis=0)
    return sig


def get_median_pairwise_correlation(E):
    C = np.corrcoef(E)
    sel = np.triu_indices(C.shape[0], k=1)
    return np.median(C[sel])


def get_signature_label(GO, sig, max_length=40):
    count = ' (%d:%d/%d)' % (sig.pc, len(sig.genes), sig.K)
    enr = sig.enr
    return GO.terms[enr.term[0]].get_pretty_format(
        omit_acc=True, max_name_length=max_length) + count


def variance_filter(genes, E, top):
    # filter genes by variance
    a = np.argsort(np.var(E, axis=1, ddof=1))[::-1]
    n = E.shape[0]
    sel = np.zeros(n, dtype=np.bool_)
    sel[a[:top]] = True
    sel = np.nonzero(sel)[0]
    genes = [genes[i] for i in sel]
    E = E[sel, :]
    return genes, E


def get_sig_matrix_figure(
        sig_matrix, max_name_length=50, include_id=False,
        highlight_sig=None, highlight_source=None,
        emin=-3.0, emax=3.0,
        font_size=12, title_font_size=16,
        margin_left=150, margin_bottom=50,
        show_sample_labels=False,
        matrix_kw=None, **kwargs):

    # colorbar label
    colorbar_label = kwargs.pop('colorbar_label', None)
    if colorbar_label is None:
        colorbar_label = 'Expression'

    heatmap = sig_matrix.get_heatmap(
        max_name_length=max_name_length,
        include_id=include_id,
        highlight_sig=highlight_sig,
        highlight_source=highlight_source,
        colorbar_label=colorbar_label,
        matrix_kw=matrix_kw,
    )

    fig = heatmap.get_figure(
        yaxis_label='Signatures',
        emin=emin, emax=emax,
        show_sample_labels=show_sample_labels,
        font_size=font_size, title_font_size=title_font_size,
        margin_left=margin_left, margin_bottom=margin_bottom,
        **kwargs
    )

    return fig

def get_sig_figure(
        sig, sig_matrix=None,
        include_pval=True,
        emin=None, emax=None,
        margin_left=70, margin_bottom=50, margin_top=50,
        show_sample_labels=False,
        #matrix_kw=None,
        sig_matrix_kw=None,
        **kwargs):

    #if matrix_kw is None:
    #    matrix_kw = {}

    if sig_matrix_kw is None:
        sig_matrix_kw = {}

    #assert isinstance(matrix_kw, dict)
    assert isinstance(sig_matrix_kw, dict)

    # generate heatmap
    heatmap = sig.get_heatmap(sig_matrix=sig_matrix,
                              sig_matrix_kw=sig_matrix_kw)

    title = sig.get_label(include_pval=include_pval)

    yaxis_label = kwargs.pop('yaxis_label', 'Genes')

    # generate figure
    fig = heatmap.get_figure(
        title=title, yaxis_label=yaxis_label,
        emin=emin, emax=emax,
        show_sample_labels=show_sample_labels,
        margin_left=margin_left,
        margin_bottom=margin_bottom, margin_top=margin_top,
        **kwargs
    )

    return fig

def read_gopca_result(path):
    """Read GO-PCA result from pickle."""
    with io.open(path, 'rb') as fh:
        G = pickle.load(fh)
    if isinstance(G, gopca.GOPCARun):
        G = G.sig_matrix
    assert isinstance(G, gopca.GOPCASignatureMatrix)
    return G


def read_go_annotations(fn):
    ann = {}
    with io.open(fn, 'rb') as fh:
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


def cluster_signatures(S, metric='correlation', method='average',
                       reverse=False):
    # hierarchical clustering of signatures
    order_rows = cluster_rows(S, metric, method, reverse)
    return order_rows


def cluster_samples(S, metric='correlation', method='average',
                    reverse=False):
    order_cols = cluster_rows(S.T, metric, method, reverse)
    return order_cols


def get_qvalues(pvals, pi_zero=1.0):
    # implements storey-tibshirani procedure for calculating q-values
    n = pvals.size
    qvals = np.empty(n, dtype=np.float64)

    # start with largest p-value
    a = np.argsort(pvals, kind='mergesort')  # stable sort
    a = a[::-1]

    s = 1
    q = 1.0
    for i in a:
        q = min(((pi_zero*pvals[i])*n)/s, q)
        qvals[i] = q
        s += 1

    return qvals
