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
from builtins import int

from copy import deepcopy

# import pytest

# from genometools.basic import GeneSetCollection
from genometools.expression import ExpMatrix
from gopca import GOPCAConfig, GOPCA


def test_basic(my_gopca):
    assert isinstance(my_gopca, GOPCA)
    assert isinstance(repr(my_gopca), str)
    assert isinstance(str(my_gopca), str)
    assert isinstance(text(my_gopca), text)
    assert isinstance(my_gopca.hash, text)

    # test members
    assert isinstance(my_gopca.configs, list)
    assert len(my_gopca.configs) > 0
    for config in my_gopca.configs:
        assert isinstance(config, GOPCAConfig)
    assert isinstance(my_gopca.matrix, ExpMatrix)

    assert isinstance(my_gopca.num_components, int)
    assert isinstance(my_gopca.pc_seed, int)
    assert isinstance(my_gopca.pc_num_permutations, int)
    assert isinstance(my_gopca.pc_zscore_thresh, float)
    assert isinstance(my_gopca.pc_max_components, int)

    # test copying
    other = deepcopy(my_gopca)
    assert other is not my_gopca
    assert other == my_gopca
    other.configs = 2*other.configs
    assert other != my_gopca

def test_simple_setup(my_gopca):
    config = my_gopca.configs[0]
    other = GOPCA.simple_setup(my_gopca.matrix,
                               config.user_params, config.gene_sets,
                               config.gene_ontology)