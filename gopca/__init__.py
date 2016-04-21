from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import pkg_resources

from gopca.config import GOPCAConfig
from gopca.signature import GOPCASignature
from gopca.result import GOPCAResult
from gopca.run import GOPCARun
from gopca.plotter import GOPCAPlotter
from gopca.go_pca import GOPCA

__version__ = str(pkg_resources.require('gopca')[0].version)

__all__ = ['GOPCAConfig', 'GOPCA', 'GOPCARun', 'GOPCAPlotter',
           'GOPCAResult', 'GOPCASignature']
