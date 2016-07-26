Preparing Gene Ontology gene sets
=================================

Gene set file format
--------------------

GO-PCA relies on gene sets derived from GO annotations (associations of genes
with specific GO terms). Therefore, one of the inputs to GO-PCA is a "GO annotation
file", which is a tab-delimited text file that contains a list of of GO terms
and the sets of genes annotated with them. An example snippet of a valid GO
annotation file is shown below:

::
    
    GO:0000795      GO      CC      synaptonemal complex    BLM,CCNB1IP1,FKBP6,MEI4,RNF212,STAG3,SYCE1,SYCE2,SYCE3,SYCP2,TEX11,UBE2I
    GO:0000796      GO      CC      condensin complex       NCAPD2,NCAPD3,NCAPG,NCAPH,SMC2,SMC4
    GO:0000808      GO      CC      origin recognition complex      LRWD1,ORC1,ORC2,ORC3,ORC4,ORC5,ORC6,REPIN1
    GO:0000812      GO      CC      Swr1 complex    ANP32E,BRD8,EP400,ING3,KAT5,RUVBL1,RUVBL2,TRRAP
    GO:0000813      GO      CC      ESCRT I complex MVB12A,MVB12B,TSG101,UBAP1,VPS28,VPS37A,VPS37B,VPS37C,VPS37D
    ...

Each row corresponds to one GO term, where the first column contains the unique
GO identifier. The second column is currently ignored by GO-PCA, and the third
column contains the "domain" of the GO term. Finally, the fourth column
contains the name of the GO term, and the last column contains a
comma-separated list of genes that are annotated with that term.