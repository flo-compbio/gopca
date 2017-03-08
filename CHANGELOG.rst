..
    Copyright (c) 2017 Florian Wagner
    
    This file is part of GO-PCA.
    
    GO-PCA is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License, Version 3,
    as published by the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.
    
    You should have received a copy of the GNU Affero General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.


Changelog
=========

Version 0.2.3 (2017-03-08)
--------------------------

- Made the following command-line scripts work again (and added tests):
  
  - gopca_print_info.py
  - gopca_extract_signatures.py
  - gopca_extract_signatures_excel.py
  - gopca_extract_signature_matrix.py
  - gopca_plot_signature_matrix.py

- Updated section on the command-line interace in the documentation, and added
  it to the table of contents.

- Internals:
  
  - Changed naming scheme of some source files (main.py => cli/go_pca.py,
    go_pca.py => gopca.py)
