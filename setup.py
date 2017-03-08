# Copyright (c) 2015-2017 Florian Wagner
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

import sys
import os
import io

from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))
description = ('GO-PCA: An Unsupervised Method to Explore Gene Expression '
               'Data Using Prior Knowledge')

# get long description from file
with io.open(path.join(here, 'README.rst'), encoding='UTF-8') as fh:
    long_description = fh.read()

install_requires = [
    'six>=1.5.2, <2',
    'future>= 0.16, <1',
    'unicodecsv>= 0.14.1, <1',
    'xlsxwriter>= 0.7.7, <1',
    'genometools>=0.2.6, <0.3',
    'setuptools>=27.2.0',
]

if sys.version_info < (3, 0):
    # We're running Python 2.x
    # => install the Python 3 configparser backport
    install_requires.append(
        'configparser>=3.2, <4',
    )

# do not require installation if built by ReadTheDocs
# (we mock these modules in docs/source/conf.py)
if 'READTHEDOCS' not in os.environ or \
        os.environ['READTHEDOCS'] != 'True':
    install_requires.extend([
        #'six>= 1.10.0, <2',
        'numpy>=1.8, <2',
        'pandas>=0.18, <1',
        'scipy>=0.14, <1',
        'scikit-learn>=0.14, <1',
        'plotly>=1.9.6, <3',
    ])
else:
    install_requires.extend([
        #'six>=1.5.2, <2',
    ])
    

setup(
    name='gopca',

    version='0.2.3',

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
        'Programming Language :: Python :: 3.5',
    ],

    keywords='unsupervised analysis gene expression data ' + \
             'transcriptomics prior knowledge',

    # packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    # packages= ['gopca', 'gopca.scripts', 'gopca.plotting'],
    packages=find_packages(exclude=['docs', 'tests*']),

    # libraries = [],

    install_requires=install_requires,

    extras_require={
        'docs': [
            'sphinx',
            'sphinx-bootstrap-theme',
            'sphinx-argparse',
            'mock'
        ],
        'tests': [
            'pytest >=2.9.1, < 4',
            'pytest-cov >=2.2.1, < 3',
            'requests >=2.10.0, <3',
        ],
    },

    # tests_require=[]

    # data
    # package_data={},

    # data outside package
    # data_files=[],

    # executable scripts
    entry_points={
        'console_scripts': [
            # GO-PCA main script
            'go-pca.py = gopca.cli.go_pca:main',

            # processing scripts
            'gopca_print_info.py = '
                'gopca.cli.print_info:main',
            'gopca_extract_signatures.py = '
                'gopca.cli.extract_signatures:main',
            'gopca_extract_signatures_excel.py = '
                'gopca.cli.extract_signatures_excel:main',
            'gopca_extract_signature_matrix.py = '
                'gopca.cli.extract_signature_matrix:main',

            #'gopca_convert_to_matlab.py = '
            #    'gopca.cli.convert_to_matlab:main',
            #'gopca_filter_signatures.py = '
            #    'gopca.cli.filter_signatures:main',
            #'gopca_combine_signatures.py = '
            #    'gopca.cli.combine_signatures:main',

            # plotting scripts
            'gopca_plot_signature_matrix.py = '
                'gopca.cli.plot_signature_matrix:main',

        ],
    },
)
