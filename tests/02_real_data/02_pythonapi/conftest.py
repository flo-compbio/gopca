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

import logging

import pytest


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@pytest.fixture(scope='session')
def my_gopca_run(my_gopca):
    logger.info('Starting GO-PCA test run...')
    gopca_run = my_gopca.run()
    return gopca_run


#@pytest.fixture(scope='session')
#def my_gopca_file(my_gopca_run, my_output_pypath):
#    gopca_file = text(my_output_pypath.join('gopca_run.pickle'))
#    my_gopca_run.write_pickle(gopca_file)
#    return gopca_file

@pytest.fixture(scope='session')
def my_sig_matrix(my_gopca_run):
    return my_gopca_run.sig_matrix
