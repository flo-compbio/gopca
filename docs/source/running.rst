Running GO-PCA
==============

This section explains the input file formats and shows the usage of the `go-pca.py` command line script.

Demo
----

Check out a `demonstration of GO-PCA <http://nbviewer.ipython.org/github/flo-compbio/gopca/blob/master/notebooks/GO-PCA_Demo.ipynb>`_ based on the `DMAP dataset <http://dx.doi.org/10.1016/j.cell.2011.01.004>`_.

Expression file format
----------------------

GO-PCA expects the expression file to be a tab-delimited text file that contains the gene expression values in a matrix layout. The first row contains the sample names, the first column represents gene names (the content of the top left cell is ignored). A mini-example of a valid expression file with only five genes and three samples is shown below:

::

    ignored Sample1 Sample2 Sample3
    IGBP1   8.64947 8.01958 7.95444
    MYC 7.61296 7.38281 7.58559
    SMAD1   8.84338 8.41662 8.94365
    MDM1    6.17908 6.07470 5.59411
    CD44    7.64093 7.56293 7.58277


GO annotation file format
-------------------------

GO-PCA relies on GO annotations, consisting of associations of genes with specific GO terms. The input file format is a tab-delimited text file, formatted as follows:

::
    
    GO:0000795      GO      CC      synaptonemal complex    BLM,CCNB1IP1,FKBP6,MEI4,RNF212,STAG3,SYCE1,SYCE2,SYCE3,SYCP2,TEX11,UBE2I
    GO:0000796      GO      CC      condensin complex       NCAPD2,NCAPD3,NCAPG,NCAPH,SMC2,SMC4
    GO:0000808      GO      CC      origin recognition complex      LRWD1,ORC1,ORC2,ORC3,ORC4,ORC5,ORC6,REPIN1
    GO:0000812      GO      CC      Swr1 complex    ANP32E,BRD8,EP400,ING3,KAT5,RUVBL1,RUVBL2,TRRAP
    GO:0000813      GO      CC      ESCRT I complex MVB12A,MVB12B,TSG101,UBAP1,VPS28,VPS37A,VPS37B,VPS37C,VPS37D

Each row corresponds to one GO term, where the first column contains the unique GO identifier. The second column is currently ignored by GO-PCA, and the third column contains the "domain" of the GO temr. Finally, the fourth column contains the name of the GO term, and the last column contains a comma-separated list of genes that are annotated with that term.


Running GO-PCA from the command line
------------------------------------

.. ".. code-block:: bash
    
    go-pca.py -g [gene_file] -a [annotation_file] -t [ontology_file] -e [expression_file] -o [output_file]

.. argparse::
   :ref: gopca.go_pca_objects.GOPCAArgumentParser
   :prog: go-pca.py
