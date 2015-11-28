import pkg_resources

from gopca.input import GOPCAInput
from gopca.signature import GOPCASignature
from gopca.go_pca import GOPCA
from gopca.output import GOPCAOutput

__version__ = pkg_resources.require('gopca')[0].version

__all__ = ['GOPCAInput','GOPCA','GOPCAOutput']
#__all__ = ['GOPCAInput','GOPCA','GOPCAOutput','GOPCASignature']
