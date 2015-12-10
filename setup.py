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

import sys
import os

from setuptools import setup, find_packages, Extension
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
description = 'GO-PCA: An Unsupervised Method to Explore Gene Expression ' + \
        'Data Using Prior Knowledge'

long_description = ''
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='gopca',

    version='1.1rc9',

    description=description,
    long_description=long_description,

    url='https://github.com/flo-compbio/gopca',

    author='Florian Wagner',
    author_email='florian.wagner@duke.edu',

    license='GPLv3',

    # see https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        'Programming Language :: Python :: 2.7',
    ],

    keywords='unsupervised analysis gene expression data transcriptomics ' + \
        'prior knowledge',

    #packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    packages= ['gopca', 'gopca.scripts', 'gopca.plotting'],
    #packages = find_packages(exclude = ['docs']),

    #libraries = [],

    install_requires=['setuptools', 'networkx', 'xlsxwriter',
            'numpy', 'scipy', 'matplotlib', 'cython', 'scikit-learn',
            'genometools>=1.2rc2', 'goparser>=1.1', 'xlmhg>=1.1rc3'],

    extras_require={
            'docs': ['sphinx','sphinx-bootstrap-theme','sphinx-argparse','mock']
    },

    # data
    #package_data={},

    # data outside package
    #data_files=[],

    # executable scripts
    entry_points={
        'console_scripts': [
            # pre-processing scripts
            'gopca_extract_go_annotations.py = gopca.extract_go_annotations:main',

            # GO-PCA scripts
            'go-pca.py = gopca.main:main',

            # processing scripts
            'gopca_extract_signatures.py = gopca.scripts.extract_signatures:main',
            'gopca_extract_signatures_excel.py = gopca.scripts.extract_signatures_excel:main',
            'gopca_extract_signature_matrix.py = gopca.scripts.extract_signature_matrix:main',

            # plotting scripts
            'gopca_plot_signature_matrix.py = gopca.plotting.plot_signature_matrix:main',
            'gopca_plot_signature.py = gopca.plotting.plot_signature:main',
        ],
    },
)
