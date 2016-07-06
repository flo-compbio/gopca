# Copyright (c) 2016 Florian Wagner
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

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import str as text
# from builtins import int as newint

from copy import deepcopy

import pytest
import numpy as np
from plotly import graph_objs as go

from genometools.basic import GeneSet
from genometools.expression import ExpMatrix
from genometools.enrichment import GSEResult
from genometools.expression.visualize import ExpHeatmap

from xlmhg import get_xlmhg_test_result
from gopca import GOPCASignature


@pytest.fixture
def my_ranked_genes():
    genes = ['a', 'c', 'd', 'f', 's']
    return genes

@pytest.fixture
def my_gene_set(my_ranked_genes):
    gene_set = GeneSet('TestID', 'TestName', my_ranked_genes)
    return gene_set

@pytest.fixture
def my_samples():
    samples = ['s1', 's2', 's3']
    return samples

@pytest.fixture
def my_v():
    v = np.uint8([1, 0, 1, 1, 0, 1] + [0] * 12 + [1, 0])
    return v

@pytest.fixture
def my_mhg_result(my_v):
    indices = np.uint16(np.nonzero(my_v)[0])
    N = my_v.size
    X = 1
    L = N
    mhg_result = get_xlmhg_test_result(N, indices, X, L)
    return mhg_result

@pytest.fixture
def my_sig_genes(my_ranked_genes, my_mhg_result):
    return my_ranked_genes[:my_mhg_result.k]

@pytest.fixture
def my_pc():
    return 1

@pytest.fixture
def my_gse_result(my_gene_set, my_mhg_result, my_ranked_genes):
    res = my_mhg_result
    gse_result = GSEResult(my_gene_set, res.N, res.indices, my_ranked_genes,
                           res.X, res.L, res.stat, res.cutoff, res.pval)
    return gse_result

@pytest.fixture
def my_matrix(my_sig_genes, my_samples):
    p = len(my_sig_genes)
    n = len(my_samples)
    X = np.arange(p*n).reshape(p, n)
    matrix = ExpMatrix(genes=my_sig_genes, samples=my_samples, X=X)
    return matrix

@pytest.fixture
def my_signature(my_pc, my_gse_result, my_matrix):
    signature = GOPCASignature(my_pc, my_gse_result, my_matrix)
    return signature

def test_basic(my_signature):
    assert isinstance(my_signature, GOPCASignature)
    assert isinstance(repr(my_signature), str)
    assert isinstance(str(my_signature), str)
    assert isinstance(text(my_signature), text)
    assert isinstance(my_signature.hash, text)

    other = deepcopy(my_signature)
    assert other is not my_signature
    assert other == my_signature
    other._pc = abs(my_signature._pc) + 1
    assert other != my_signature

def test_label(my_signature):
    assert isinstance(my_signature.label, text)
    label = my_signature.get_label(max_name_length=5, include_stats=True,
                                   include_id=True, include_pval=True,
                                   include_coll=True)
    assert isinstance(label, text)

def test_heatmap(my_gopca_signature):
    heatmap = my_gopca_signature.get_heatmap(
        sig_matrix_kw={'cluster_samples': False},
        colorbar_label=r'Median-centered expression (log<sub>2</sub>-RPKM)',
    )
    assert isinstance(heatmap, ExpHeatmap)

    fig = heatmap.get_figure(
        emin=-3, emax=3,
        height=800, font_size=12, title_font_size=18, show_sample_labels=True,
        margin_bottom=100,
    )
    assert isinstance(fig, go.graph_objs.Figure)