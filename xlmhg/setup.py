from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
#from Cython.Distutils import Extension

import numpy as np

include_dirs = []
library_dirs = []

ext_modules = []
ext_modules.append(Extension("xlmHG_cython", ["xlmHG_cython.pyx"],\
		library_dirs=library_dirs,
		extra_link_args=['-fPIC'],
		include_dirs=[np.get_include()]+include_dirs))

setup(
  name = 'cython tools',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules,
)
