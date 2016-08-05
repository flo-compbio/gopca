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
from builtins import str as text

from copy import deepcopy

import pytest

from genometools.basic import GeneSetCollection
from gopca import GOPCAParams, GOPCAConfig

import six
if six.PY2:
    import cPickle as pickle
else:
    import pickle


def test_basic(my_config):
    assert isinstance(my_config, GOPCAConfig)
    assert isinstance(repr(my_config), str)
    assert isinstance(str(my_config), str)
    assert isinstance(text(my_config), text)
    assert isinstance(my_config.hash, text)

    # test members
    assert isinstance(my_config.user_params, GOPCAParams)
    assert isinstance(my_config.gene_sets, GeneSetCollection)
    assert isinstance(my_config.gene_ontology, type(None))
    assert isinstance(my_config.params, type(None))

    # test copying
    other = deepcopy(my_config)
    assert other is not my_config
    assert other == my_config
    other.final_params = other.params
    assert other != my_config