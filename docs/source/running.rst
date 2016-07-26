Running *GO-PCA*
================

This section describes the expression file format used by GO-PCA, and documents
the ``go-pca.py`` command that is used to run GO-PCA. For information on how
to provide the gene sets to GO-PCA, please see the previous section,
`gene_sets`.

The expression file format
--------------------------

The main input to GO-PCA is an expression matrix, with rows representing genes,
and columns representing samples. GO-PCA expects the expression matrix to be described in a tab-delimited text file that contains the gene expression values in a matrix layout. The first row contains the sample names, and the first column represents gene names (the content of the top left cell is ignored). A mini-example of a valid expression file with only five genes and three samples is shown below:

::

    ignored Sample1 Sample2 Sample3
    IGBP1   8.64947 8.01958 7.95444
    MYC     7.61296 7.38281 7.58559
    SMAD1   8.84338 8.41662 8.94365
    MDM1    6.17908 6.07470 5.59411
    CD44    7.64093 7.56293 7.58277


The GO-PCA workflow
-------------------

Currently, all GO-PCA functionalities are accessible through the command-line.
This allows users that are not familiar with the Python programming language to
make full use of GO-PCA. Since many different analyses and visualizations can
be performed based on the results of a single GO-PCA run, the GO-PCA workflow
is organized as follows:

1. GO-PCA is run using the command ``go-pca.py`` All results are
   written to a binary *result file*.

2. A variety of `command-line scripts <cli>` can be used to extract detailed
   information about the signatures generated in different formats. These
   scripts use the result file as *input* data. In this way, GO-PCA does not
   need to be re-run for every operation.

.. _go_pca:

Running GO-PCA: ``go-pca.py``
-----------------------------

.. ".. code-block:: bash
    
    go-pca.py -g [gene_file] -a [annotation_file] -t [ontology_file] -e [expression_file] -o [output_file]

``go-pca.py`` is the command to run *GO-PCA*. All parameters can either be
spcefied directly on the command line, or in a separate configuration file,
using the ``-c`` option.

.. note::

  The configuration file is expected to follow the Windows "INI-style" format,
  with a single "[GO-PCA]" section, followed by "parameter=value" entries. 
  If a configuration file is given, and a parameter is set both in the
  configuration file and on the command line, the command line setting takes
  precedence.

The only required parameters are:::

 -e  (The expression file.)
 -s  (The gene set file.)
 -o  (The output file.)

However, if the expression matrix is not pre-filtered to only contain expressed
genes, it is also highly advisable to specify the ``-G`` option.

.. argparse::
   :ref: gopca.main.get_argument_parser
   :prog: go-pca.py
