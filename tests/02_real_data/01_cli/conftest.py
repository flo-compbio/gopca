# Copyright (c) 2017 Florian Wagner
#
# This file is part of GO-PCA.
#
# GO-PCA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
# from builtins import *
_oldstr = str
from builtins import open
from builtins import str

import logging
import subprocess as subproc

import pytest
import six

if six.PY2:
    import cPickle as pickle
else:
    import pickle


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@pytest.fixture(scope='session')
def my_gopca_file(
        my_expression_file,
        my_fly_gene_set_file,
        my_gene_ontology_file,
        my_output_pypath):
    """Run the go-pca.py script to generate a pickle file with the results."""

    gopca_file = str(my_output_pypath.join('gopca.pickle'))

    p = subproc.Popen(
        'go-pca.py -o %s -e %s -s %s -t %s -G 8000'
        % (gopca_file, my_expression_file, my_fly_gene_set_file,
           my_gene_ontology_file),
        shell=True, stdout=subproc.PIPE, stderr=subproc.PIPE)

    stdout, stderr = p.communicate()
    if p.returncode != 0:
        print('Stderr:')
        for l in stderr.decode('utf-8').split('\n'):
            print(l)
    assert p.returncode == 0, str(stderr)  # no errors

    return gopca_file


@pytest.fixture(scope='session')
def my_gopca_run(my_gopca_file):
    """Load the results (a `GOPCARun` object) from the pickle file."""
    with open(my_gopca_file, 'rb') as fh:
        gopca_run = pickle.load(fh)
    return gopca_run


#@pytest.fixture(scope='session')
#def my_gopca_file(my_gopca_run, my_output_pypath):
#    gopca_file = text(my_output_pypath.join('gopca_run.pickle'))
#    my_gopca_run.write_pickle(gopca_file)
#    return gopca_file
