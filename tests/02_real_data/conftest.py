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
import shutil

import pytest
import requests
import logging

from genometools.expression import ExpMatrix
from genometools.expression.filter import filter_variance
from genometools.basic import GeneSetCollection
from genometools.ontology import GeneOntology

from gopca import GOPCAParams, GOPCAConfig, GOPCA

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def download_file(url, path):
    """Downloads a file."""
    r = requests.get(url, stream=True)
    if r.status_code == 200:
        with open(path, 'wb') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)


@pytest.fixture(scope='session')
def my_params():
    params = GOPCAParams()
    # params.set_param('sel_var_genes', 8000)
    params.set_param('mHG_X_frac', 0.25)
    params.set_param('mHG_X_min', 5)
    params.set_param('mHG_L', 1000)
    return params


@pytest.fixture(scope='session')
def my_data_pypath(tmpdir_factory):
    pypath = tmpdir_factory.mktemp('gopca_data', numbered=False)
    return pypath


@pytest.fixture(scope='session')
def my_output_pypath(tmpdir_factory):
    pypath = tmpdir_factory.mktemp('gopca_output', numbered=False)
    return pypath


@pytest.fixture(scope='session')
def my_expression_file(my_data_pypath):
    logger.info('Starting download of expression file...')
    url = r'https://www.dropbox.com/s/6lqqdp9meeyasin/fly_timecourse_expression.tsv?dl=1'
    path = text(my_data_pypath.join('fly_timecourse_expression.tsv'))
    download_file(url, path)
    return path


@pytest.fixture(scope='session')
def my_gene_ontology_file(my_data_pypath):
    logger.info('Starting download of gene ontology file...')
    url = r'https://www.dropbox.com/s/rvox0wi96162it2/go-basic_2016-01-18.obo?dl=1'
    path = text(my_data_pypath.join('go-basic_2016-01-18.obo'))
    download_file(url, path)
    return path


@pytest.fixture(scope='session')
def my_fly_gene_set_file(my_data_pypath):
    """Drosophila gene set file."""
    logger.info('Starting download of fly gene set file...')
    url = r'https://www.dropbox.com/s/rcu4qqvilrg68ko/GO_gene_sets_fly_ensembl83_goa54_ontology2016-01-18.tsv?dl=1'
    path = text(my_data_pypath.join('GO_gene_sets_fly_ensembl83_goa54_ontology2016-01-18.tsv'))
    download_file(url, path)
    return path


@pytest.fixture(scope='session')
def my_config(my_params, my_gene_ontology_file, my_fly_gene_set_file):
    gene_ontology = GeneOntology.read_obo(my_gene_ontology_file)
    gene_sets = GeneSetCollection.read_tsv(my_fly_gene_set_file)
    config = GOPCAConfig(my_params, gene_sets, gene_ontology)
    return config


@pytest.fixture(scope='session')
def my_matrix_filtered(my_expression_file):
    matrix = ExpMatrix.read_tsv(my_expression_file)
    matrix_filtered = filter_variance(matrix, 8000)
    return matrix_filtered


@pytest.fixture(scope='session')
def my_gopca(my_config, my_matrix_filtered):
    configs = [my_config]
    gopca = GOPCA(my_matrix_filtered, configs, pc_seed=123456789)
    return gopca
