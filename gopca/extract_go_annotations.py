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
GO annotation data from the `UniProt-GOA database`__, extracting a list of all
genes annotated with each GO term, and stores the results in a tab-delimited
text file.

The five columns of the output file are:

1. GO term ID (e.g., "GO\:0006260").
2. "GO" (constant).
3. Abbreviated domain of the GO term (e.g., "BP" for
   ``biological_process``).
4. GO term name (e.g., "DNA replication").
5. Comma-separated list of genes (represented by their gene symbols) associated
   with the GO term.


__ uniprot_goa_

.. _uniprot_goa: http://www.ebi.ac.uk/GOA

Examples
--------

-   Extract the GO annotations from UniProt-GOA release 149, for all
    human protein coding genes from `Ensembl`__ release 82, retaining only GO terms
    that have at least 5 and no more than 200 genes annotated with them.

    __ ensembl_

    In a first step, extract a list of all protein-coding genes from the
    `Ensembl GTF file`__, using the script `extract_protein_coding_genes.py` from
    the :mod:`genometools` package:


    __ gtf_file

    .. code-block:: bash

        $ extract_protein_coding_genes.py \\
            -a Homo_sapiens.GRCh38.82.gtf.gz \\
            -o protein_coding_genes_human.tsv

    In the second step, extract the GO annotations, based on the human
    `gene association file `__ (in GAF 2.0 format) from UniProt-GOA, and the
    corresponding version of the `gene ontology file`__ (in OBO 1.2 format)
    from the Gene Ontology Consortium:

    __ gaf_file_

    .. code-block:: bash

        $ gopca_extract_go_annotations.py -g protein_coding_genes_human.tsv \\
            -t go-basic.obo -a gene_association.goa_human.149.gz \\
            --min-genes-per-term 5 --max-genes-per-term 200 \\
            -o go_annotations_human.tsv

.. _ensembl: http://www.ensembl.org
.. _gtf_file: ftp://ftp.ensembl.org/pub/release-82/gtf/homo_sapiens/Homo_sapiens.GRCh38.82.gtf.gz
.. _gaf_file: ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.149.gz
.. _obo_file: http://viewvc.geneontology.org/viewvc/GO-SVN/ontology-releases/2015-10-12/go-basic.obo?revision=29122

"""

# we don't assume that gene names are sorted

import sys
import os

import argparse
import csv
import gzip
import logging
import cPickle as pickle

import numpy as np
import networkx as nx

from genometools import misc
from goparser import GOParser
#from gopca import common

def get_argument_parser():
    """Function to obtain the argument parser.

    Returns
    -------
    A fully configured `argparse.ArgumentParser` object.

    Notes
    -----
    This function is used by the `sphinx-argparse` extension for sphinx.
    """
    
    description = """Script for determining the set of genes annotated with \
            each GO term."""

    parser = argparse.ArgumentParser(description = description)

    # input files
    parser.add_argument('-g','--gene-file',required=True,
            help="""Path of tab-delimited file containing all "valid" gene
                    symbols.""")

    parser.add_argument('-t','--ontology-file',required=True,
            help='Path of ontology file (in OBO format).')

    parser.add_argument('-a','--gene-association-file',required=True,
            help='Path of gene association file (in GAF format).')

    # output file
    parser.add_argument('-o','--output-file',required=True,
            help='Path of output file.')

    # evidence
    parser.add_argument('-e','--select-evidence',nargs='+',default=[],
            help="""List of three-letter evidence codes to include.
                    If not specified, include all evidence types.""")

    # which GO terms to include in final output?
    parser.add_argument('--min-genes-per-term',type=int,default=0,
            help="""Exclude GO terms that have fewer than the specified number
                    of genes annotated with them. Disabled (0) by default.""")

    parser.add_argument('--max-genes-per-term',type=int,default=0,
            help="""Exclude GO terms that have more than the specified number
                    of genes annotated with them. Disabled (0) by default.""")

    # logging options
    parser.add_argument('-l','--log-file',default=None,
            help='Path of log file.')

    parser.add_argument('-q','--quiet',action='store_true',
            help='Suppress all output except warnings and errors.')

    parser.add_argument('-v','--verbose',action='store_true',
            help='Enable verbose output. Ignored if ``--quiet`` is specified.')

    # legacy options
    parser.add_argument('--part-of-cc-only',action='store_true',
            help="""If enabled, ignore ``part_of`` relations outside the
                    ``cellular_component`` (CC) domain.""")

    return parser

def main(args=None):
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
    ontology_file = args.ontology_file
    gene_association_file = args.gene_association_file
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
    log_level = logging.INFO
    if quiet:
        log_level = logging.WARNING
    elif args.verbose:
        log_level = logging.DEBUG

    logger = misc.configure_logger('', log_file = log_file,
            log_level = log_level)

    # checks
    assert os.path.isfile(gene_file)
    assert os.path.isfile(ontology_file)
    assert os.path.isfile(gene_association_file)

    # read genes and sort them
    genes = sorted(misc.read_single(args.gene_file))
    n = len(genes)
    logger.info('Read %d genes.', n); sys.stdout.flush()

    # Read GO term definitions and parse UniProtKB GO annotations
    if len(select_evidence) == 1 and (not select_evidence[0].strip(' ')):
        select_evidence = []

    misc.configure_logger('goparser', log_file = log_file,
            log_level = log_level)
    GO = GOParser()
    GO.parse_ontology(ontology_file,part_of_cc_only=False)
    GO.parse_annotations(gene_association_file,gene_file,select_evidence=select_evidence)

    #with open(go_pickle_file) as fh:
    #   GO = pickle.load(fh)

    # Get sorted list of GO term IDs
    all_term_ids = sorted(GO.terms.keys())

    logger.info('Obtaining GO term associations...')
    n = len(all_term_ids)
    term_gene_counts = []
    term_ids = []
    term_genes = []
    for j,id_ in enumerate(all_term_ids):
        tg = GO.get_goterm_genes(id_)
        assert isinstance(tg,set)
        c = len(tg)
        if c >= min_genes and c <= max_genes:
            term_gene_counts.append(c)
            term_ids.append(id_)
            term_genes.append(tg)
    term_gene_counts = np.int64(term_gene_counts)

    # remove GO terms that are perfectly redundant, keep descendant terms
    logger.info('Testing for perfect overlap...')
    m = len(term_ids)
    #genesets = [set(np.nonzero(A[:,j])[0]) for j in range(m)]
    #term_gene_count = np.sum(A,axis=0,dtype=np.int64)
    G = nx.Graph()
    G.add_nodes_from(range(m))
    for j1 in range(m):
        #if (j1+1) % 1000 == 0: print j1+1, ; sys.stdout.flush()
        c = term_gene_counts[j1]
        tg = term_genes[j1]
        for j2 in range(m):
            if j2 >= j1: break
            if c == term_gene_counts[j2] and tg == term_genes[j2]:
                G.add_edge(j1,j2)

    sel = np.ones(m,dtype=np.bool_)
    affected = 0
    for k,cc in enumerate(nx.connected_components(G)):
        if len(cc) == 1: # singleton
            continue
        affected += len(cc)
        for j1 in cc:
            keep = True
            term = GO.terms[term_ids[j1]]
            for j2 in cc:
                if j1 == j2: continue
                if term_ids[j2] in term.descendants:
                    keep = False
                    break
            if not keep:
                sel[j1] = False
    logger.info('# affected terms: %d', affected)
    logger.info('# perfectly redundant descendant terms: %d', np.sum(np.invert(sel)))

    sel = np.nonzero(sel)[0]
    term_ids = [term_ids[j] for j in sel]
    term_genes = [term_genes[j] for j in sel]
    logger.info('Selected %d / %d non-redundant GO terms.', sel.size,m)

    # write output file
    logger.info('Writing output file...')
    p = len(genes)
    with open(output_file,'w') as ofh:
        writer = csv.writer(ofh,dialect='excel-tab',lineterminator='\n',quoting=csv.QUOTE_NONE)
        for j,(id_,tg) in enumerate(zip(term_ids,term_genes)):
            term = list(GO.terms[id_].get_tuple())
            writer.writerow(term + [','.join(sorted(tg))])
    logger.info('done!');

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
