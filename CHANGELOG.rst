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


Version 1.2 (work in progress)
------------------------------

- Full Python 3.5 support
- Integrated code with genometools 2.0
- Cleaned up code
- Full Python API
- Uses plotly as plotting backend (dropped matplotlib)
- Reorganized internal folder structure (combined all command-line scripts in
  one folder)

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
