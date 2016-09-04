..
    Copyright (c) 2015, 2016 Florian Wagner
    
    This file is part of GO-PCA.
    
    GO-PCA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License, Version 3,
    as published by the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.


Changelog
=========

Version 2.1.0 (2016-08-28)
--------------------------

- Improved algorithm for generating signatures based on significantly enriched
  gene sets, and introduced a new parameter, ``sig_min_genes`` to specify the
  minimum number of genes in each signature. Effectively, this makes gene set
  enrichment analysis and signature generation more independent of each other,
  and results in more coherent signatures (i.e., genes in a signature should
  have more similar expression patterns with the new algorithm).

- Improved the "Installing GO-PCA" section of the documentation. The section
  now includes accurate instructions for Windows and Mac OS X.

- Reduced the verbosity of GO-PCA output. The old verbosity can be restored
  by passing ``verbose=True`` to GO-PCA, or by running the ``go-pca.py`` script
  with the ``-v`` option.

- Changed data structure used for representing signature matrices
  (`GOPCASignatureMatrix` now inherits from `genometools.expression.ExpMatrix`).

- Added simplified functions for plotting heatmaps of signature matrices
  and signature genes to the Python API: `GOPCASignatureMatrix.get_figure()`
  and `GOPCASignature.get_figure()`.

2.1.1 (2016-09-04)
~~~~~~~~~~~~~~~~~~

- Travis deploy step now generates and uploads wheels.
- Fixed bug in `GOPCASignatureMatrix.filter_collection_signatures()`.

2.1.2 (2016-09-04)
~~~~~~~~~~~~~~~~~~

- Changed genometools dependency to version 2.0.2.

Version 2.0.0 (2016-08-01)
--------------------------

- Generalized the API to support multiple "configurations" as input. Each
  configuration (`GOPCAConfig`) consists of a set of GO-PCA parameters
  (`GOPCAParams`), a list of gene sets
  (`genometools.basic.GeneSetCollection`), and, optionally, a gene ontology
  (`genometools.ontology.GeneOntology`). GO-PCA is performed independently
  for each configuration, and the signatures generated are merged into a
  single signature matrix. The old API is available as `GOPCA.init_simple()`.

Version 2.0.1 (2016-08-16)
~~~~~~~~~~~~~~~~~~~~~~~~~~
- Some bugfixes

Version 1.2.0 (2016-07-26)
--------------------------

- Full Python 3.5 support (on top of Python 2.7)
- New Python API for running GO-PCA, accessing results and generating figures
- Integrated code with genometools 2.0
- Cleaned up code
- Added tests (using py.test)
- Switched to plotly as plotting backend (dropped matplotlib)

1.2.1 (2016-07-26)
~~~~~~~~~~~~~~~~~~

- Updated changelog, fix version number in documentation

1.2.4 (2016-08-04)
~~~~~~~~~~~~~~~~~~

- Fixed bug that made it impossible to run GO-PCA from the command line
  (`go-pca.py`)

1.2.5 (2016-08-05)
~~~~~~~~~~~~~~~~~~

- Fixed a bug that prevented the Gene Ontology file from being loaded when
  running GO-PCA from the command line (`go-pca.py`)

Version 1.1.3 (2016-02-23)
--------------------------

- Added a lot more documentation (see http://gopca.readthedocs.org).
- Added ``gopca_plot_all_signatures.py`` command.

Version 1.1.2 (2016-02-04)
--------------------------

- Fixed a bug in the local filter that resulted in incorrect filtering

Version 1.1.0 (2016-02-02)
--------------------------
- Restructured the documentation (work still ongoing).
- The "GO annotation file" is now more accurately being referred to as
  "gene set file"
- Modified file format for supplying gene sets (added a "description" column)
- GO-PCA is intended to work with arbitrary gene sets, and the code now
  reflects this as well (parameter names, variable names, etc.)
- Added support for INI-style GO-PCA configuration file ("go-pca.py -c ...")
- Improved documentation of command-line pamaraters ("go-pca.py --help")
- Improved support for config and result MD5 hash values
  (they are now shown by "gopca_print_info.py")

Version 1.1rc12
---------------
- Reorganized internal class structure (see `GOPCARun` in ``run.py``)
- Added Unicode support
- Added ``gopca_print_info.py`` script to quickly inspect GO-PCA output data

Version 1.1rc10
---------------

- Fixed a bug that occurred when no ontology file is given.

Version 1.1rc9
--------------

"Visible" changes:

- Changed the names of some command line arguments.
- Added/improved documentation of command line parameters for all scripts.
- Changed the sphinx html documentation theme to bootstrap
  ("sphinx-bootstrap-theme").

Internal changes:

- Improved code documentation.
- The GO-PCA parameters `disable_local_filter` and `disable_global_filter` are
  now called `no_nocal_filter` and `no_global_filter`, respectively.
- The expression matrix is now represented using the `ExpMatrix` class from the
  `genometools` package.
- Shared parameter for plotting scripts are now obtained using functions from
  the `plotting.params` module (this greatly reduced code redundancy).
- Loggers are no longer class members, and are instead defined as global
  variables within each module. This is consistent with the recommended naming
  scheme that uses `logging.getLogger(__name__)`, thus naming a logger after
  the module. This helped to simplify the class structures.
