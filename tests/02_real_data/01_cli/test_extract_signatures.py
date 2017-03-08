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

"""Tests for the `gopca_extract_signatures.py` script."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
# from builtins import *
from builtins import open
from builtins import str as text

import subprocess as subproc

import pytest


@pytest.fixture(scope='module')
def my_output_file(my_output_pypath):
    """The output file."""
    return text(my_output_pypath.join('gopca_signatures.tsv'))


def test_no_error(my_gopca_file, my_output_file):
    """Run the script and make sure that the return code is zero."""
    p = subproc.Popen(
        'gopca_extract_signatures.py -g %s -o %s'
        % (my_gopca_file, my_output_file),
        shell=True, stdout=subproc.PIPE, stderr=subproc.PIPE)

    stdout, stderr = p.communicate()
    if p.returncode != 0:
        print('Stderr:')
        for l in stderr.decode('utf-8').split('\n'):
            print(l)
    assert p.returncode == 0, str(stderr)  # no errors
