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

root = 'gopca'

here = path.abspath(path.dirname(__file__))
description = 'GO-PCA: An Unsupervised Method to Explore Biological Heterogeneity in Gene Expression Data Using Prior Knowledge'

long_description = ''
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='gopca',

    version='1.1rc1',

    description=description,
    long_description=long_description,

    url='https://github.com/flo-compbio/gopca',

    author='Florian Wagner',
    author_email='florian.wagner@duke.edu',

    license='GPLv3',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        'Programming Language :: Python :: 2.7',
    ],

    keywords='unsupervised analysis gene expression data transcriptomics prior knowledge',

    #packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    packages=[root],

	#libraries = [],

    install_requires=['cython','numpy','scipy','scikit-learn','networkx','genometools','goparser','xlmhg','sphinx','sphinx-rtd-theme'],

	# development dependencies
    #extras_require={},

	# data
    #package_data={},

	# data outside package
    #data_files=[],

	# executable scripts
    entry_points={
        'console_scripts': [
            'extract_go_annotations.py = gopca.extract_go_annotations:main',
            'go-pca.py = gopca.go_pca:main',
			'extract_signatures.py = gopca.scripts.extract_signatures:main',
			'extract_signature_matrix.py = gopca.scripts.extract_signature_matrix:main',
            'extract_matlab_file.py = gopca.scripts.extract_matlab_file:main',
            'plot_signature_matrix.py = gopca.plotting.plot_signature_matrix:main',
            'plot_signature_correlation_matrix.py = gopca.plotting.plot_signature_correlation_matrix:main',
			'plot_within_signature_correlations.py = gopca.plotting.plot_within_signature_correlations:main',
            'plot_term_by_pc_matrix.py = gopca.plotting.plot_term_by_pc_matrix:main',
			'gopca_test_components.py = gopca.test_components:main'
        ],
    },
)
