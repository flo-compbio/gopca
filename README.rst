GO-PCA
======

GO-PCA is an unsupervised method to **explore gene expression data using prior
knowledge**. The `GO-PCA paper`__ has recently been published in PLOS ONE.

__ go_pca_paper_

GO-PCA combines `principal component analysis (PCA)`__  with
`nonparametric GO enrichment analysis`__ in order to generate **signatures**,
i.e., small sets of genes that are both strongly correlated and closely
functionally related. It then produces a **signature matrix**, designed to
serve as a systematic and easily interpretable representation of biologically
relevant expression patterns.

__ pca_
__ go_enrich_

.. _go_pca_paper: https://dx.doi.org/10.1371/journal.pone.0143196
.. _pca: https://en.wikipedia.org/wiki/Principal_component_analysis
.. _go_enrich: https://dx.doi.org/10.1186/1471-2105-10-48

`Documentation <https://gopca.readthedocs.org/en/latest>`_
----------------------------------------------------------

- `Demo <http://nbviewer.ipython.org/github/flo-compbio/gopca/blob/master/notebooks/GO-PCA_Demo.ipynb>`_
- `Installation <https://gopca.readthedocs.org/en/latest/install.html>`_
- `Running GO-PCA <https://gopca.readthedocs.org/en/latest/running.html>`_

How to Cite GO-PCA
------------------

If you use GO-PCA in your research, please cite `Wagner (PLOS ONE, 2015)`__

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
