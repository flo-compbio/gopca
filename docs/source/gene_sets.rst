Obtaining Gene Sets
===================

GO-PCA requires two types of inputs: First, a gene expression matrix (see
`Running GO-PCA <running>`), and second, a list of gene sets.
These gene sets represent the prior knowledge used by GO-PCA, and form the
basis for signature generation. Generally speaking, they should represent
functional categories of potential relevance to the expression dataset
analyzed. This section discusses different types of gene sets and how to
obtain gene set files for use with GO-PCA.

In short, pre-generated GO-derived gene set files (and corresponding Gene
Ontology files in ``OBO`` format) are available for `download`__ (see below).
Alternatively, users can create and use their own, customized gene sets.

__ gene_sets_
.. # __ download_gene_sets_
.. # __ create_gene_sets_

.. #contents:: Contents
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

Up-to-date, GO-derived gene sets based on high-confidence GO evidence codes are
`available for download`__. For more details on the procedure used to derive
these gene sets from GO annotation data, please see the included README.

__ gene_sets_

Note that GO-PCA should be provided with both the gene set file (.tsv) and the
matching Gene Ontology file (.obo), in order to be able to apply all filters
(see previous :ref:`Note <filtering_note>`). For more information on how to
specify both gene sets and Gene Ontology when running GO-PCA from the command-
line, see :doc:`running`.

.. _gene_sets: https://www.dropbox.com/sh/ovekr0h7l60onoa/AACNMWUQOJnxdatLge205IFUa?dl=0

Generating custom GO-derived gene sets
--------------------------------------

Custom GO-derived gene sets can be generated using the script
:ref:`gopca_extract_go_gene_sets.py <extract_go_gene_sets>`.