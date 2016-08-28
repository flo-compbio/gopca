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
# from builtins import *
from builtins import open
from builtins import str as text

# import os
# import shutil

import pytest
import logging

import numpy as np

from xlmhg import get_xlmhg_test_result
from genometools.expression import ExpMatrix, ExpProfile
from genometools.basic import GeneSet, GeneSetCollection
from genometools.enrichment import RankBasedGSEResult
# from genometools.ontology import GeneOntology

from gopca import GOPCAParams, GOPCAConfig, \
                  GOPCASignature, GOPCASignatureMatrix, \
                  GOPCA

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)


@pytest.fixture(scope='session')
def my_params():
    params = GOPCAParams()
    return params


@pytest.fixture(scope='session')
def my_gene_sets():
    gene_sets = GeneSetCollection([
        GeneSet('GeneSet1', 'First gene set', ['a', 'b', 'd'],
                 source='TestSource', collection='TestCollection',
                 description='The first test GeneSet.'),
        GeneSet('GeneSet2', 'Second gene set', ['a', 'c', 'd'],
                source='TestSource', collection='TestCollection',
                description='The second test GeneSet.'),
    ])
    return gene_sets


@pytest.fixture(scope='session')
def my_config(my_params, my_gene_sets):
    config = GOPCAConfig(my_params, my_gene_sets)
    return config


@pytest.fixture(scope='session')
def my_matrix():
    genes = ['a', 'b', 'c', 'd', 'e', 'f']
    samples = ['s1', 's2', 's3']
    X = np.arange(18, dtype=np.float64).reshape(6, 3)
    matrix = ExpMatrix(genes=genes, samples=samples, X=X)
    return matrix


@pytest.fixture(scope='session')
def my_v():
    v = np.uint8([1,1,1,0,0,0])
    return v


@pytest.fixture(scope='session')
def my_gopca(my_config, my_matrix):
    configs = [my_config]
    my_gopca = GOPCA(my_matrix, configs)
    return my_gopca


@pytest.fixture(scope='session')
def my_rank_based_result(my_matrix, my_v):

    indices = np.uint16(np.nonzero(my_v)[0])
    ind_genes = [my_matrix.genes[i] for i in indices]
    gs_genes = list(ind_genes)
    gene_set = GeneSet(genes=gs_genes, id='Random1', name='Random gene Set 1')
    N = my_v.size
    X = 1
    L = N
    ## stat, n_star, pval = xlmhg_test(my_v, X, L)
    res = get_xlmhg_test_result(N, indices, X, L)
    result = RankBasedGSEResult(gene_set, N, indices, ind_genes, X, L,
                                res.stat, res.cutoff, res.pval)
    return result


# TO-DO: GOPCASignature, GOPCASignatureMatrix, GOPACRun
@pytest.fixture(scope='session')
def my_signature(my_matrix, my_rank_based_result):
    sig_genes = my_rank_based_result.genes_above_cutoff
    seed = ExpProfile(my_matrix.loc[sig_genes].mean(axis=0))
    sig = GOPCASignature(1, my_rank_based_result, seed,
                         my_matrix.loc[sig_genes])
    return sig


@pytest.fixture(scope='session')
def my_other_signature(my_matrix, my_rank_based_result):
    sig_genes = my_rank_based_result.genes_above_cutoff[:-1]
    seed = ExpProfile(my_matrix.loc[sig_genes].mean(axis=0))
    sig = GOPCASignature(0, my_rank_based_result, seed,
                         my_matrix.loc[sig_genes])
    return sig


@pytest.fixture(scope='session')
def my_sig_matrix(my_signature, my_other_signature):
    signatures = [my_signature, my_other_signature]
    sig_matrix = GOPCASignatureMatrix.from_signatures(signatures)
    return sig_matrix