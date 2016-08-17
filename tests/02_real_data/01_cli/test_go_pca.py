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

import pytest

import os
import subprocess as subproc

@pytest.fixture(scope='session')
def my_output_file(my_output_pypath):
    return text(my_output_pypath.join('command-line_output.pickle'))


def is_writable(path):
    try:
        with open(path, 'a'):
            pass
    except:
        return False
    return True


def test_script(my_expression_file,
                my_fly_gene_set_file,
                my_gene_ontology_file,
                my_output_file):

    assert is_writable(my_output_file)

    p = subproc.Popen(
        'go-pca.py -o %s -e %s -s %s -t %s -D 5'
        % (my_output_file, my_expression_file, my_fly_gene_set_file,
           my_gene_ontology_file),
        shell=True, stdout=subproc.PIPE, stderr=subproc.PIPE)

    stdout, stderr = p.communicate()
    print('Stderr:')
    for l in stderr.decode('utf-8').split('\n'):
        print(l)
    assert p.returncode == 0, str(stderr)  # no errors