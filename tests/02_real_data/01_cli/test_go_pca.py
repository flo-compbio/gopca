# Copyright (c) 2016, 2017 Florian Wagner
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

"""Tests for the main script, `go-pca.py`."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
# from builtins import *
from builtins import open
from builtins import str as text

import os
import subprocess as subproc

import pytest

from gopca import GOPCARun, GOPCASignatureMatrix


#@pytest.fixture(scope='module')
#def my_output_file(my_output_pypath):
#    return text(my_output_pypath.join('command-line_output.pickle'))


def test_correct_output(my_gopca_run):
    """Test if the GO-PCA main script produced the correct output."""
    assert isinstance(my_gopca_run, GOPCARun)

    sig_matrix = my_gopca_run.sig_matrix
    assert isinstance(sig_matrix, GOPCASignatureMatrix)
    assert len(sig_matrix.signatures) == 51
