from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import pkg_resources

from .config import GOPCAConfig
from .signature import GOPCASignature
from .signature_matrix import GOPCASignatureMatrix
from .run import GOPCARun
from .go_pca import GOPCA

__version__ = str(pkg_resources.require('gopca')[0].version)

__all__ = ['GOPCAConfig', 'GOPCA', 'GOPCARun',
           'GOPCASignatureMatrix', 'GOPCASignature']
