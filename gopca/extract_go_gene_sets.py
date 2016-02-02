#!/usr/bin/env python2.7

# Copyright (c) 2015 Florian Wagner
#
# This file is part of GO-PCA.
#
# GO-PCA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Script for determining the set of genes annotated with each GO term.

This script (see `main` function) uses the :mod:`goparser` package to parse
GO annotation data from the `UniProt-GOA database`__, extracts a list of all
genes annotated with each GO term, and stores the results in a tab-delimited
text file.

The output file format is specified in the
`genometools.basic.geneset.GeneSetDB` class, and contains six colums:

1. Gene set ID - here: GO term ID (e.g., "GO\:0006260").
2  Source - here: Gene Ontology ("GO").
3. Collection - here: GO domain (e.g., "BP" for ``biological_process``).
4. Gene set name - here: GO term name (e.g., "DNA replication").
5. Comma-separated list of genes - here: all genes annotated with the GO term.
6. Gene set description - here: definition of the GO term

__ uniprot_goa_

.. _uniprot_goa: http://www.ebi.ac.uk/GOA

Examples
--------

-   Extract gene set basde on the GO annotations from UniProt-GOA release 149,
    for all human protein coding genes from `Ensembl`__ release 82, retaining
    only GO terms that have at least 5 and no more than 200 genes annotated
    with them.

    __ ensembl_

    In a first step, extract a list of all protein-coding genes from the
    `Ensembl GTF file`__, using the script `extract_protein_coding_genes.py`
    from the :mod:`genometools` package:

    __ gtf_file

    .. code-block:: bash

        $ extract_protein_coding_genes.py \\
            -a Homo_sapiens.GRCh38.82.gtf.gz \\
            -o protein_coding_genes_human.tsv

    In the second step, extract the gene sets, based on the human
    `GO annotation file `__ (in GAF 2.0 format) from UniProt-GOA, and the
    corresponding version of the `gene ontology file`__ (in OBO 1.2 format)
    from the Gene Ontology Consortium:

    __ gaf_file_

    .. code-block:: bash

        $ gopca_extract_go_gene_sets.py -g protein_coding_genes_human.tsv \\
            -t go-basic.obo -a gene_association.goa_human.149.gz \\
            --min-genes-per-term 5 --max-genes-per-term 200 \\
            -o go_gene_sets_human.tsv

.. _ensembl: http://www.ensembl.org
.. _gtf_file: ftp://ftp.ensembl.org/pub/release-82/gtf/homo_sapiens/Homo_sapiens.GRCh38.82.gtf.gz
.. _gaf_file: ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.149.gz
.. _obo_file: http://viewvc.geneontology.org/viewvc/GO-SVN/ontology-releases/2015-10-12/go-basic.obo?revision=29122

"""

import sys
import os
import argparse
import gzip
import logging
import cPickle as pickle
import textwrap

import numpy as np
import networkx as nx

import unicodecsv as csv

from genometools import misc
from goparser import GOParser
from gopca import util
from gopca import cli

def get_argument_parser():
    """Function to obtain the argument parser.

    Returns
    -------
    A fully configured `argparse.ArgumentParser` object.

    Notes
    -----
    This function is used by the `sphinx-argparse` extension for sphinx.
    """
    #prog = 'gopca_extract_go_annotations.py'
    desc = """Write a file containing the genes annotated with each GO term."""

    parser = cli.get_argument_parser(desc = desc)

    file_mv = cli.file_mv
    str_type = cli.str_type

    # input files
    g = parser.add_argument_group('Input and output files')

    g.add_argument('-g', '--gene-file', required = True,
            type = str_type, metavar = file_mv,
            help="""Path of tab-delimited file containing all "valid" gene
                    symbols.""")

    g.add_argument('-t', '--gene-ontology-file', required = True,
            metavar = cli.file_mv,
            help='Path of ontology file (in OBO format).')

    g.add_argument('-a', '--go-annotation-file', required = True,
            metavar = cli.file_mv,
            help='Path of GO annotation file (in GAF format).')

    # output file
    g.add_argument('-o', '--output-file', required = True,
            metavar = cli.file_mv,
            help='Path of output file.')

    g = parser.add_argument_group('Other options')

    # evidence
    g.add_argument('-e', '--select-evidence', nargs = '*', default = None,
            metavar = '<evidence code, ...>',
            help="""List of three-letter evidence codes to include.
                    If not specified, include all evidence types.""")

    # which GO terms to include in final output?
    g.add_argument('--min-genes-per-term', type = int, default=0,
            metavar = cli.int_mv,
            help="""Exclude GO terms that have fewer than the specified number
                    of genes annotated with them. Disabled (0) by default.""")

    g.add_argument('--max-genes-per-term', type = int, default=0,
            metavar = cli.int_mv,
            help="""Exclude GO terms that have more than the specified number
                    of genes annotated with them. Disabled (0) by default.""")

    # legacy options
    g.add_argument('--part-of-cc-only', action = 'store_true',
            help="""If enabled, ignore ``part_of`` relations outside the
                    ``cellular_component`` (CC) domain.""")

    # reporting options
    cli.add_reporting_args(parser)

    return parser

def main(args = None):
    """Extract GO annotations and store in tab-delimited text file.

    Parameters
    ----------
    args: argparse.Namespace object, optional
        The argument values. If not specified, the values will be obtained by
        parsing the command line arguments using the `argparse` module.

    Returns
    -------
    int
        Exit code (0 if no error occurred).
 
    """


    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    gene_file = args.gene_file
    gene_ontology_file = args.gene_ontology_file
    go_annotation_file = args.go_annotation_file
    output_file = args.output_file

    select_evidence = args.select_evidence
    min_genes = args.min_genes_per_term
    max_genes = args.max_genes_per_term

    part_of_cc_only = args.part_of_cc_only

    # logging parameters
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # configure root logger
    logger = util.get_logger(log_file = log_file, quiet = quiet,
            verbose = verbose)

    # checks
    assert os.path.isfile(gene_file)
    assert os.path.isfile(gene_ontology_file)
    assert os.path.isfile(go_annotation_file)

    # read genes and sort them
    genes = sorted(misc.read_single(args.gene_file))
    n = len(genes)
    logger.info('Read %d genes.', n)

    # Read GO term definitions and parse UniProtKB GO annotations
    GO = GOParser()

    logger.info('Parsing gene ontology file...')
    GO.parse_ontology(gene_ontology_file, part_of_cc_only = False)

    logger.info('Parsing GO annotation file...')
    genes = misc.read_single(gene_file)
    GO.parse_annotations(go_annotation_file, genes,
            select_evidence = select_evidence)

    D = GO.get_gene_sets(min_genes = min_genes, max_genes = max_genes)
    D.write_tsv(output_file)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
