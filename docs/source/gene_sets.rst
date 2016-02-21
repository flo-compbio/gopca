Obtaining Gene Sets
===================

Besides a gene expression matrix (see `Running GO-PCA <running>`_), GO-PCA
requires a list of gene sets. These gene sets represent the prior knowledge
used by GO-PCA, and form the basis for signature generation. Generally
speaking, they should represent functional categories of potential relevance
to the expression dataset analyzed. This section discusses different types of
gene sets and how to obtain gene set files for use with GO-PCA.

Pre-generated GO-derived gene set files (and corresponding Gene Ontology files
in ``OBO`` format) are available for download (`see below`__). Users can also
`create`__ and use their own, customized gene sets.

__ download_gene_sets_
__ create_gene_sets_

.. contents:: Contents
    :depth: 2
    :local:
    :backlinks: none


GO-derived vs. custom gene sets
-------------------------------

GO-PCA was originally designed to work with gene sets derived from Gene
Ontology (GO) annotation data in the `UniProt-GOA`__ database. Each *GO term*
gives rise to a particular gene set, which is defined as the set of all genes
*annotated* with that term. However, not all such gene sets are necessarily
useful for GO-PCA: There are many overly specific GO terms that are associated
with very few genes, and there are some overly general GO terms that are
associated with too many genes. Moreover, many GO annotations are uncertain,
sometimes based solely on computational predictions. These issues need to be
taken into consideration when deriving gene sets from GO annotation data.

__ uniprot_goa_

.. _uniprot_goa: http://www.ebi.ac.uk/GOA

.. _filtering_note:

.. note::
    
    The hierarchical structure of the Gene Ontology also results in many
    highly overlapping gene sets. However, this is anticipated by the *GO-PCA*
    algorithm, which includes two different filtering steps designed to limit
    the generation of redundant signatures. For this filtering to work as
    intended, GO-PCA also needs to be given an "ontology file", containing the
    structure of the Gene Ontology (in ``OBO`` format).


While GO-derived gene sets represent a high-quality, all-purpose repository
of prior knowledge, many GO-PCA analyses might benefit from alternative or
additional gene sets that either serve to supplement the UniProt-GOA database,
or represent more context-specific prior knowledge that a researcher might have
acquired from previous experiments. GO-PCA will treat such *custom* gene sets
in the same way as GO-derived gene sets, with the minor exception that certain
certain filtering rules do not apply (see previous :ref:`Note
<filtering_note>`).

.. _go_pca_paper: https://dx.doi.org/10.1371/journal.pone.0143196

.. _download_gene_sets:

Downloading pre-generated gene sets
-----------------------------------

Up-to-date, high-quality GO-derived gene sets for multiple species are
`available for download`__. For more details on the procedure used to derive
these gene sets from GO annotation data, please see the included README.

__ gene_sets_

Note that GO-PCA should be provided with both the gene set file (.tsv) and the
matching Gene Ontology file (.obo), in order to be able to apply all filters
(see previous :ref:`Note <filtering_note>`). For more information on how to
specify both gene sets and Gene Ontology when running GO-PCA from the command-
line, see :doc:`running`.

.. _gene_sets: https://www.dropbox.com/sh/m0r7uqnfdr5x0xu/AADqqJ-8VzPchBRhDm50QxWaa?dl=0

.. _create_gene_sets:

Creating gene sets from scratch
-------------------------------

Depending on their goals, users have two options for creating their own gene
sets. If they seek to derive gene sets from GO annotations using their own
parameters (e.g., making their own selection of which GO evidence codes to
include), they can use the script `gopca_extract_go_gene_sets.py` (see below).
In contrast, if they wish to generate gene sets independent of GO annotations,
they have to create their own gene set files.

The gene set file format
~~~~~~~~~~~~~~~~~~~~~~~~

GO-PCA expects gene sets to be defined in a tab-delimited text file,
with each row corresponding to one gene set. There are six columns, defined as
follows:

1. **Gene Set ID** - A unique identifier for the gene set. In the case of
   GO-based gene sets, this should be a GO term ID.
2. **Source** - A string indicating the source of the gene set (e.g., "GO" for
   Gene Ontology).
3. **Collection** - A string indicating the sub-category or group the gene set
   belongs to (e.g., "CC" for the "cellular component" domain of the Gene
   Ontology).
4. **Name** - The name of the gene set.
5. **Genes** - A comma-separated list of genes that are in the gene set.
6. **Description** - A short description of the gene set.

Currently, GO-PCA ignores fields 2 and 6, but it may use the information in
those fields in future releases. Here is a small example of a valid gene set
file containing only five gene sets:

::
    
    GO:0000795      GO      CC      synaptonemal complex    BLM,CCNB1IP1,FKBP6,MEI4,RNF212,STAG3,SYCE1,SYCE2,SYCE3,SYCP2,TEX11,UBE2I    A proteinaceous scaffold found between homologous chromosomes during meiosis.
    GO:0000796      GO      CC      condensin complex       NCAPD2,NCAPD3,NCAPG,NCAPH,SMC2,SMC4 A multisubunit protein complex that plays a central role in chromosome condensation.
    GO:0000808      GO      CC      origin recognition complex      LRWD1,ORC1,ORC2,ORC3,ORC4,ORC5,ORC6,REPIN1  A multisubunit complex that is located at the replication origins of a chromosome.
    GO:0000812      GO      CC      Swr1 complex    ANP32E,BRD8,EP400,ING3,KAT5,RUVBL1,RUVBL2,TRRAP A multisubunit protein complex that is involved in chromatin remodeling. It is required for the incorporation of the histone variant H2AZ into chromatin. In S. cerevisiae, the complex contains Swr1p, a Swi2/Snf2-related ATPase, and 12 additional subunits.
    GO:0000813      GO      CC      ESCRT I complex MVB12A,MVB12B,TSG101,UBAP1,VPS28,VPS37A,VPS37B,VPS37C,VPS37D    An endosomal sorting complex required for transport. It consists of the class E vacuolar protein sorting (Vps) proteins and interacts with ubiquitinated cargoes.


Deriving gene sets from GO annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: gopca.extract_go_gene_sets.get_argument_parser
   :prog: gopca_extract_go_gene_sets.py
