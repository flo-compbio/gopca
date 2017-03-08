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


Running GO-PCA from the command-line
------------------------------------

For users who are unfamiliar with Python, the most convenient way of running GO-PCA
is through the command line. In general, this can be done in a two-step process:

1. GO-PCA is run using the command `go-pca.py <cli>` All results are
   written to a binary *result file*.

2. A variety of additional `command-line scripts <cli>` can be used to extract
   detailed information about the signatures generated, and to plot the GO-PCA
   signature matrix.
   These scripts use the result file as *input* data. In this way, GO-PCA does
   not need to be re-run for every operation.