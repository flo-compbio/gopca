from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import pkg_resources

from .config import GOPCAConfig
from .signature import GOPCASignature
from .result import GOPCAResult
from .run import GOPCARun
from .plotter import GOPCAPlotter
from .go_pca import GOPCA

__version__ = str(pkg_resources.require('gopca')[0].version)

__all__ = ['GOPCAConfig', 'GOPCA', 'GOPCARun', 'GOPCAPlotter',
           'GOPCAResult', 'GOPCASignature']
