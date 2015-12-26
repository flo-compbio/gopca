# Copyright (c) 2015 Florian Wagner
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

import pkg_resources

from gopca.config import GOPCAConfig
from gopca.signature import GOPCASignature
from gopca.output import GOPCAOutput
from gopca.run import GOPCARun
from gopca.go_pca import GOPCA

__version__ = pkg_resources.require('gopca')[0].version

__all__ = ['GOPCAInput','GOPCA','GOPCAOutput', 'GOPCARun']
