# Copyright (c) 2015, 2016 Florian Wagner
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

from __future__ import print_function

# import sys
import os
import io

import six

from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))
description = ('GO-PCA: An Unsupervised Method to Explore Gene Expression '
               'Data Using Prior Knowledge')

# get long description from file
with io.open(path.join(here, 'README.rst'), encoding='UTF-8') as fh:
    long_description = fh.read()

setup_requires = [
    'six >= 1.10.0, < 2',
]

install_requires = [
    'future >= 0.15.2, < 1',
    'six >= 1.10.0, < 2',
    'unicodecsv >= 0.14.1, < 1',
    'xlsxwriter >= 0.7.7, < 1',
]

if six.PY2:
    install_requires.append(
        'configparser >= 3.2, < 4',
    )

# do not require installation if built by ReadTheDocs
# (we mock these modules in docs/source/conf.py)
if 'READTHEDOCS' not in os.environ or \
        os.environ['READTHEDOCS'] != 'True':
    install_requires.extend([
        'numpy >= 1.8, < 2',
        'pandas >= 0.18, < 1',
        'scipy >= 0.14, < 1',
        'scikit-learn >= 0.14, < 1',
        'matplotlib >= 1.4.3, < 2',
        'plotly >= 1.9.6, < 2',
        'genometools >= 2.0rc1',
        'goparser >= 1.2rc1, < 2',
        'xlmhg >= 2.0.6, < 3',
    ])

setup(
    name='gopca',

    version='2.0rc1',

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
        'Programming Language :: Python :: 3',
    ],

    keywords='unsupervised analysis gene expression data ' + \
             'transcriptomics prior knowledge',

    # packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    # packages= ['gopca', 'gopca.scripts', 'gopca.plotting'],
    packages=find_packages(exclude=['docs', 'tests*']),

    # libraries = [],

    setup_requires=setup_requires,

    install_requires=install_requires,

    extras_require={
         'docs': ['sphinx', 'sphinx-bootstrap-theme', 'sphinx-argparse',
                  'mock']
    },

    # data
    # package_data={},

    # data outside package
    # data_files=[],

    # executable scripts
    entry_points={
        'console_scripts': [
            # pre-processing scripts
            'gopca_extract_go_gene_sets.py = '
                'gopca.extract_go_gene_sets:main',

            # GO-PCA main script
            'go-pca.py = gopca.main:main',

            # processing scripts
            'gopca_extract_signatures.py = '
                'gopca.scripts.extract_signatures:main',

            'gopca_extract_signatures_excel.py = '
                'gopca.scripts.extract_signatures_excel:main',

            'gopca_extract_signature_matrix.py = '
                'gopca.scripts.extract_signature_matrix:main',

            'gopca_convert_to_matlab.py = '
                'gopca.scripts.convert_to_matlab:main',
            'gopca_filter_signatures.py = '
                'gopca.scripts.filter_signatures:main',
            'gopca_combine_signatures.py = '
                'gopca.scripts.combine_signatures:main',
            'gopca_print_info.py = '
                'gopca.scripts.print_info:main',

            # plotting scripts
            'gopca_plot_signature_matrix.py = '
                'gopca.plotting.plot_signature_matrix:main',

            'gopca_plot_signature.py = '
                'gopca.plotting.plot_signature:main',

            'gopca_plot_all_signatures.py = '
                'gopca.plotting.plot_all_signatures:main',

            'gopca_plot_term_by_pc_matrix.py = '
                'gopca.plotting.plot_term_by_pc_matrix:main',

            # 'gopca_plot_pc_scores.py = gopca.plotting.plot_pc_scores:main',
        ],
    },
)
