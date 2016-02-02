Obtaining a Gene Set File
=========================

Besides a gene expression matrix (see `Running GO-PCA<running>`_), GO-PCA requires a list of gene sets, typically defined based on Gene Ontology (GO) annotations. Up-to-date, high-quality GO-based gene set files for multiple species are `available for download`__. The rest of this section describes the gene set file format used by GO-PCA. It is only important to understand this format if you want to use GO-PCA with a customized list of gene sets.

__ gene_sets_

.. _gene_sets: https://www.dropbox.com/sh/m0r7uqnfdr5x0xu/AADqqJ-8VzPchBRhDm50QxWaa?dl=0

The gene set file format
------------------------

GO-PCA expects gene sets to be defined in a tab-delimited plain-text file, with each row corresponding to one gene set. There are six columns, defined as follows:

1. **Gene Set ID** - A unique identifier for the gene set. In the case of GO-based gene sets, this should be a GO term ID.
2. **Source** - A string indicating the source of the gene set (e.g., "GO" for Gene Ontology).
3. **Collection** - A string indicating the sub-category or group the gene set belongs to (e.g., "CC" for the "cellular component" domain of the Gene Ontology).
4. **Name** - The name of the gene set.
5. **Genes** - A comma-separated list of genes that are in the gene set.
6. **Description** - A short description of the gene set.

Currently, GO-PCA ignores fields 2 and 6, but it may use the information in those fields in future releases. Here is a small example of a valid gene set file containing only five gene sets:

::
    
    GO:0000795      GO      CC      synaptonemal complex    BLM,CCNB1IP1,FKBP6,MEI4,RNF212,STAG3,SYCE1,SYCE2,SYCE3,SYCP2,TEX11,UBE2I    A proteinaceous scaffold found between homologous chromosomes during meiosis.
    GO:0000796      GO      CC      condensin complex       NCAPD2,NCAPD3,NCAPG,NCAPH,SMC2,SMC4 A multisubunit protein complex that plays a central role in chromosome condensation.
    GO:0000808      GO      CC      origin recognition complex      LRWD1,ORC1,ORC2,ORC3,ORC4,ORC5,ORC6,REPIN1  A multisubunit complex that is located at the replication origins of a chromosome.
    GO:0000812      GO      CC      Swr1 complex    ANP32E,BRD8,EP400,ING3,KAT5,RUVBL1,RUVBL2,TRRAP A multisubunit protein complex that is involved in chromatin remodeling. It is required for the incorporation of the histone variant H2AZ into chromatin. In S. cerevisiae, the complex contains Swr1p, a Swi2/Snf2-related ATPase, and 12 additional subunits.
    GO:0000813      GO      CC      ESCRT I complex MVB12A,MVB12B,TSG101,UBAP1,VPS28,VPS37A,VPS37B,VPS37C,VPS37D    An endosomal sorting complex required for transport. It consists of the class E vacuolar protein sorting (Vps) proteins and interacts with ubiquitinated cargoes.

