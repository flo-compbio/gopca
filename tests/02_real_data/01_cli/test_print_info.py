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
from builtins import open
from builtins import str as text

import pytest

import os
import subprocess as subproc

@pytest.fixture(scope='session')
def my_gopca_file(my_gopca_run, my_output_pypath):
    gopca_file = text(my_output_pypath.join('gopca_run.pickle'))
    my_gopca_run.write_pickle(gopca_file)
    return gopca_file


def is_writable(path):
    try:
        with open(path, 'a'):
            pass
    except:
        return False
    return True


def test_script(my_gopca_file):

    p = subproc.Popen(
        'gopca_print_info.py -g %s'
        % (my_gopca_file),
        shell=True, stdout=subproc.PIPE, stderr=subproc.PIPE)

    stdout, stderr = p.communicate()
    print('Stderr:')
    for l in stderr.decode('utf-8').split('\n'):
        print(l)
    assert p.returncode == 0, str(stderr)  # no errors