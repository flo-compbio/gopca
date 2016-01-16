..
    Copyright (c) 2015 Florian Wagner
    
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

GO-PCA
======

|docs-latest| |docs-develop|

GO-PCA (`Wagner, 2015`__) is an unsupervised method to **explore gene
expression data using prior knowledge**. This is a free and open-source
implementation of GO-PCA in Python.

__ go_pca_paper_

Briefly, GO-PCA combines `principal component analysis (PCA)`__  with
`nonparametric GO enrichment analysis`__ in order to generate **signatures**,
i.e., small sets of genes that are both strongly correlated and closely
functionally related. It then visualizes the expression profiles of all
signatures in a **signature matrix**, designed to serve as a systematic and
easily interpretable representation of biologically relevant expression
patterns.

__ pca_
__ go_enrich_

.. _go_pca_paper: https://dx.doi.org/10.1371/journal.pone.0143196
.. _pca: https://en.wikipedia.org/wiki/Principal_component_analysis
.. _go_enrich: https://dx.doi.org/10.1186/1471-2105-10-48

Documentation
-------------

- `Homepage <https://gopca.readthedocs.org/en/latest>`_
- `"DMAP" Demo <http://nbviewer.ipython.org/url/gopca.readthedocs.org/en/latest/_downloads/Demo_DMAP.ipynb>`_
- `Getting Started <https://gopca.readthedocs.org/en/latest/getting_started.html>`_

Support and Development
-----------------------

- For feature requests and bug reports, please create an `issue`__ on GitHub.
- For technical questions, please feel free to `email`__.
- If you want to contribute code to GO-PCA, please `email`__ and/or create a
  pull request on GitHub.
- For a list of the latest changes, please see the
  `changelog <changelog.rst>`_.

__ github_issue_
__ email_
__ email_

.. _github_issue: https://github.com/flo-compbio/gopca/issues
.. _email: mailto:florian.wagner@duke.edu

How to Cite GO-PCA
------------------

If you use GO-PCA in your research, please cite `Wagner (PLoS One, 2015)`__

__ wagner_pone_

.. _wagner_pone: https://dx.doi.org/10.1371/journal.pone.0143196

Copyright and License
---------------------

Copyright (c) 2015 Florian Wagner

::

  GO-PCA is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License, Version 3,
  as published by the Free Software Foundation.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

.. |docs-latest| image:: https://readthedocs.org/projects/gopca/badge/?version=latest
    :alt: Documentation Status (master branch)
    :scale: 100%
    :target: https://gopca.readthedocs.org/en/latest/?badge=latest

.. |docs-develop| image:: https://readthedocs.org/projects/gopca/badge/?version=develop
    :alt: Documentation Status (develop branch)
    :scale: 100%
    :target: https://gopca.readthedocs.org/en/develop/?badge=develop

