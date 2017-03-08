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

import os

import pytest

# import requests
# import shutil

from genometools.expression import ExpMatrix
from genometools.basic import GeneSetCollection
from genometools.ontology import GeneOntology

from gopca import GOPCARun, GOPCASignatureMatrix


def test_download(my_expression_file, my_gene_ontology_file,
                  my_fly_gene_set_file):
    """Test if required data files were downloaded successfully."""

    # expression file
    print(my_expression_file)
    assert os.path.isfile(my_expression_file)
    matrix = ExpMatrix.read_tsv(my_expression_file)
    assert isinstance(matrix, ExpMatrix)
    assert matrix.hash == 'aa7cc5e6e04d34e65058f059bcdfe5ea'

    # gene ontology file
    print(my_gene_ontology_file)
    assert os.path.isfile(my_gene_ontology_file)
    # hash not stable?
    #ontology = GeneOntology.read_obo(my_gene_ontology_file)
    #assert isinstance(ontology, GeneOntology)
    #assert ontology.hash == '978546899cfb0196ac2005d4b177725f'

    # gene set file
    print(my_fly_gene_set_file)
    assert os.path.isfile(my_fly_gene_set_file)
    gene_sets = GeneSetCollection.read_tsv(my_fly_gene_set_file)
    assert isinstance(gene_sets, GeneSetCollection)
    assert gene_sets.hash == '78b4b27e9658560a8e5993154d3228fa'


def test_gopca_run_correct_results(my_gopca_run):
    """Test if the GO-PCA test run produced the correct results."""
    assert isinstance(my_gopca_run, GOPCARun)

    sig_matrix = my_gopca_run.sig_matrix
    assert isinstance(sig_matrix, GOPCASignatureMatrix)
    assert len(sig_matrix.signatures) == 51
