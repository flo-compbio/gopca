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
description = 'GO-PCA: An Unsupervised Method to Explore Gene Expression Data Using Prior Knowledge'

long_description = ''
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='gopca',

    version='1.1rc5',

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
    packages=['gopca','gopca.scripts','gopca.plotting'],

	#libraries = [],

    install_requires=['cython','numpy','scipy','scikit-learn','networkx','genometools','goparser','xlmhg','sphinx','sphinx-rtd-theme','sphinx-argparse'],

	# development dependencies
    #extras_require={},

	# data
    #package_data={},

	# data outside package
    #data_files=[],

	# executable scripts
    entry_points={
        'console_scripts': [
			'gopca_test_components.py = gopca.test_components:main',
			'gopca_extract_go_annotations.py = gopca.extract_go_annotations:main',
			'go-pca.py = gopca.go_pca:main',

			'gopca_extract_signatures.py = gopca.scripts.extract_signatures:main',
			'gopca_extract_signatures_excel.py = gopca.scripts.extract_signatures_excel:main',
			'gopca_extract_signature_matrix.py = gopca.scripts.extract_signature_matrix:main',
			'gopca_extract_matlab_file.py = gopca.scripts.extract_matlab_file:main',

			'gopca_plot_signature_matrix.py = gopca.plotting.plot_signature_matrix:main',
			'gopca_plot_signature_correlation_matrix.py = gopca.plotting.plot_signature_correlation_matrix:main',
			'gopca_plot_within_signature_correlations.py = gopca.plotting.plot_within_signature_correlations:main',
			'gopca_plot_term_by_pc_matrix.py = gopca.plotting.plot_term_by_pc_matrix:main',
			'gopca_plot_signature.py = gopca.plotting.plot_signature:main',
			'gopca_plot_pc_variance_explained.py = gopca.plotting.plot_pc_variance_explained:main',
			'gopca_plot_pc_scores.py = gopca.plotting.plot_pc_scores:main',

			'bootstrap-go-pca.py = gopca.bootstrap_go_pca:main',
			'gopca_extract_bootstrap_sample.py = gopca.scripts.extract_bootstrap_sample:main',
			'gopca_plot_bootstrap_numbers.py = gopca.plotting.plot_bootstrap_numbers:main',
			'gopca_plot_bootstrap_sample_size_summary.py = gopca.plotting.plot_bootstrap_sample_size_summary:main',
			'gopca_plot_bootstrap_sample_size_matrix.py = gopca.plotting.plot_bootstrap_sample_size_matrix:main',
			'gopca_plot_bootstrap_pc_matrix.py = gopca.plotting.plot_bootstrap_pc_analysis:main',
			'gopca_plot_bootstrap_signature_detection.py = gopca.plotting.plot_bootstrap_signature_detection:main',
        ],
    },
)
