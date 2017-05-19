#!/usr/bin/env python

# Copyright (c) 2015-2017 Florian Wagner
#
# This file is part of GO-PCA.
#
# GO-PCA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Script for determining the set of genes annotated with each GO term.

This script (see `main` function) uses the :mod:`genometools` package to parse
GO annotation data from the `UniProt-GOA database`__, extracts a list of all
genes annotated with each GO term, and stores the results in a tab-delimited
text file.

The output file format is specified in the
`genometools.basic.geneset.GeneSetCollection` class, and contains six colums:

1. Gene set ID - here: GO term ID (e.g., "GO:0006260").
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
    `Ensembl GTF file`__, using the script
    `ensembl_extract_protein_coding_genes.py` from the :mod:`genometools`
    package:

    __ gtf_file

    .. code-block:: bash

        $ ensembl_extract_protein_coding_genes.py \\
            -a Homo_sapiens.GRCh38.88.gtf.gz \\
            -o protein_coding_genes_human.tsv

    In the second step, extract the gene sets, based on the human
    `GO annotation file `__ (in GAF 2.0 format) from UniProt-GOA, and the
    corresponding version of the `gene ontology file`__ (in OBO 1.2 format)
    from the Gene Ontology Consortium:

    __ gaf_file_

    .. code-block:: bash

        $ gopca_extract_go_gene_sets.py -g protein_coding_genes_human.tsv \\
            -t go-basic.obo -a goa_human_167.gaf.gz \\
            --min-genes-per-term 5 --max-genes-per-term 200 \\
            -o go_gene_sets_human.tsv

.. _ensembl: http://www.ensembl.org
.. _gtf_file: ftp://ftp.ensembl.org/pub/release-82/gtf/homo_sapiens/Homo_sapiens.GRCh38.82.gtf.gz
.. _gaf_file: ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.149.gz
.. _obo_file: http://viewvc.geneontology.org/viewvc/GO-SVN/ontology-releases/2015-10-12/go-basic.obo?revision=29122

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import sys
import os
import textwrap
import gzip

from genometools import misc
from genometools.expression import ExpGenome
from genometools import ontology
from genometools.ontology import GeneOntology
from genometools.basic import GeneSetCollection

from . import arguments


def get_argument_parser():
    """Function to obtain the argument parser.

    Returns
    -------
    A fully configured `argparse.ArgumentParser` object.

    Notes
    -----
    This function is used by the `sphinx-argparse` extension for sphinx.
    """

    prog = 'gopca_extract_go_gene_sets.py'
    description = 'Extract GO-derived gene sets.'
    parser = arguments.get_argument_parser(prog, description)

    file_mv = arguments.file_mv
    # name_mv = arguments.name_mv
    int_mv = arguments.int_mv
    # float_mv = arguments.float_mv
    # str_mv = arguments.str_mv

    # input files
    g = parser.add_argument_group('Input and output files')

    g.add_argument(
        '-g', '--gene-file', type=str, required=True, metavar=file_mv,
        help=textwrap.dedent("""\
            File containing list of protein-coding genes (generated using
            the script ``ensembl_extract_protein_coding_genes.py``)."""))

    g.add_argument(
        '-t', '--gene-ontology-file', type=str, required=True, metavar=file_mv,
        help='Path of ontology file (in OBO format).')

    g.add_argument(
        '-a', '--goa-association-file', type=str, required=True, metavar=file_mv,
        help='Path of UniProt-GOA Gene Association file (in GAF format).')

    # output file
    g.add_argument(
        '-o', '--output-file', type=str, required=True, metavar=file_mv,
        help='Path of output file.')

    g = parser.add_argument_group('Other options')

    # evidence
    g.add_argument(
        '-e', '--evidence-codes', nargs='*',
        default=['IDA', 'IGI', 'IMP', 'ISO', 'ISS', 'IC', 'NAS', 'TAS'],
        metavar='<evidence code, ...>',
        help=textwrap.dedent("""\
            List of three-letter evidence codes to include.
            If empty, include all evidence types.
            [IDA, IGI, IMP, ISO, ISS, IC, NAS, TAS]"""))

    # which GO terms to include in final output?
    g.add_argument(
        '--min-genes-per-term', type=int, default=5, metavar=int_mv,
        help=textwrap.dedent("""\
            Exclude GO terms that have fewer than the specified number
            of genes annotated with them. Set to 0 to disable. [5]"""))

    g.add_argument(
        '--max-genes-per-term', type=int, default=200, metavar=int_mv,
        help=textwrap.dedent("""\
            Exclude GO terms that have more than the specified number
            of genes annotated with them. Set to 0 to disable. [200]"""))

    # legacy options
    g.add_argument(
        '--part-of-cc-only', action='store_true',
        help=textwrap.dedent("""\
            If enabled, ignore ``part_of`` relations outside the
            ``cellular_component`` (CC) domain."""))

    # reporting options
    arguments.add_reporting_args(parser)

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
 
    Raises
    ------
    SystemError
        If the version of the Python interpreter is not >= 2.7.
    """
    vinfo = sys.version_info
    if not (vinfo >= (2, 7)):
        raise SystemError('Python interpreter version >= 2.7 required, '
                          'found %d.%d instead.' % (vinfo.major, vinfo.minor))

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    gene_file = args.gene_file
    gene_ontology_file = args.gene_ontology_file
    goa_association_file = args.goa_association_file
    output_file = args.output_file

    evidence_codes = args.evidence_codes
    min_genes = args.min_genes_per_term
    max_genes = args.max_genes_per_term

    part_of_cc_only = args.part_of_cc_only

    # logging parameters
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # configure root logger
    logger = misc.get_logger(log_file=log_file, quiet=quiet,
                             verbose=verbose)

    logger.info('Selected evidence codes: %s', ', '.join(evidence_codes))
    logger.info('Min. number of genes per gene set: %d', min_genes)
    logger.info('Max. number of genes per gene set: %d', max_genes)

    # checks
    assert os.path.isfile(gene_file)
    assert os.path.isfile(gene_ontology_file)
    assert os.path.isfile(goa_association_file)

    # configure root logger
    log_stream = sys.stdout
    if output_file == '-':
        # if we print output to stdout, redirect log messages to stderr
        log_stream = sys.stderr

    logger = misc.get_logger(log_stream=log_stream, log_file=log_file,
                             quiet=quiet, verbose=verbose)

    # extract protein-coding genes from Ensembl GTF file
    exp_genome = ExpGenome.read_tsv(gene_file)

    # parse Gene Ontology
    gene_ontology = GeneOntology.read_obo(gene_ontology_file)

    # parse UniProt-GOA gene association file
    with gzip.open(goa_association_file, 'rt', encoding='ascii') as fh:
        go_annotations = ontology.parse_gaf(
            fh, gene_ontology, ev_codes=evidence_codes, genome=exp_genome)

    # extract GO-based gene sets
    gene_sets = ontology.get_goa_gene_sets(go_annotations)
    logger.info('Generated %d GO-derived gene sets', len(gene_sets))

    # filter gene sets based on size
    if min_genes > 0:
        old_size = len(gene_sets)
        gene_sets = GeneSetCollection(
            gs for gs in gene_sets if gs.size >= min_genes)
        logger.info('Excluded %d gene sets with too few genes.',
                    old_size-len(gene_sets))

    if max_genes > 0:
        old_size = len(gene_sets)
        gene_sets = GeneSetCollection(
            gs for gs in gene_sets if gs.size <= max_genes)
        logger.info('Excluded %d gene sets with too many genes.',
                    old_size-len(gene_sets))
        

    # writing output file
    gene_sets.write_tsv(output_file)
    logger.info('Wrote %s GO-derived gene sets to output file "%s".',
                len(gene_sets), output_file)

    return 0


if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
