Running GO-PCA
==============

This section describes the expression file format used by GO-PCA, and documents the `go-pca.py` command that is used to run GO-PCA.

Expression file format
----------------------

The main input to GO-PCA is an expression matrix, with rows representing genes, and columns representing samples. GO-PCA expects the expression matrix to be described in a tab-delimited text file that contains the gene expression values in a matrix layout. The first row contains the sample names, and the first column represents gene names (the content of the top left cell is ignored). A mini-example of a valid expression file with only five genes and three samples is shown below:

::

    ignored Sample1 Sample2 Sample3
    IGBP1   8.64947 8.01958 7.95444
    MYC     7.61296 7.38281 7.58559
    SMAD1   8.84338 8.41662 8.94365
    MDM1    6.17908 6.07470 5.59411
    CD44    7.64093 7.56293 7.58277



Running GO-PCA from the command line
------------------------------------

To run GO-PCA, use the ``go-pca.py`` command.

.. ".. code-block:: bash
    
    go-pca.py -g [gene_file] -a [annotation_file] -t [ontology_file] -e [expression_file] -o [output_file]

.. argparse::
   :ref: gopca.main.get_argument_parser
   :prog: go-pca.py
