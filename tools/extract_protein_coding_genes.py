#!/usr/bin/env python

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

"""
Extracts protein coding genes from GENCODE annotations.
"""

import sys
import os
import gzip
import csv
import re
import gzip
import argparse

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-a','--annotation-file',required=True)
	parser.add_argument('-o','--output-file',required=True)
	parser.add_argument('-e','--exclude-chromosomes',default=[],nargs='+')

	parser.add_argument('-i','--case-insensitive',action='store_true')

	return parser.parse_args()

def open_plain_or_gzip(fn):
	try:
		gzip.open(fn).next()
		return gzip.open(fn)
	except IOError:
		return open(fn)

attr_sep = re.compile(r"(?<!\\)\s*;\s*") # use negative lookbehind to make sure we don't split on escaped semicolons ("\;")
def parse_attributes(s):
	''' parses the 9th field (attributes) of a GFF/GTF entry into a dictionary '''
	attr = {}
	atts = attr_sep.split(s)
	for a in atts:
		#print a
		kv = a.split(' ')
		if len(kv) == 2:
			k,v = kv
			v = v.strip('"')
			attr[k] = v
	return attr	

def read_chromlen(fn):
	CL = {}
	with open(fn) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		for l in reader:
			CL[l[0]] = int(l[1])
	return CL

def main(args):

	exclude_chromosomes = args.exclude_chromosomes
	case_insensitive = args.case_insensitive
	print "Excluded:",exclude_chromosomes; sys.stdout.flush()
	genes = {}
	gene_status = {}
	transcripts = {}
	sources = {}
	chromosomes = {}
	i = 0
	n_t = 0
	n_g = 0
	with gzip.open(args.annotation_file) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		for l in reader:
			if len(l) == 1: continue # skip header
			i += 1
			if i % int(1e6) == 0: print i # report progress

			if l[2] == 'gene':
				attr = parse_attributes(l[8])
				if attr['gene_type'] in ['protein_coding','polymorphic_pseudogene']:
					chrom = l[0]
					source = l[1]
					id_ = attr['gene_id']
					name = attr['gene_name']
					if case_insensitive: name = name.upper()
					status = attr['gene_status']

					# convert to Ensembl
					if chrom.startswith('chr'): chrom = chrom[3:]
					if chrom == 'M': chrom = 'MT'

					chromosomes[id_] = chrom
					genes[id_] = name
					sources[id_] = source
					gene_status[id_] = status
					n_g += 1

			elif l[2] == 'transcript':
				attr = parse_attributes(l[8])
				if attr['transcript_type'] == 'protein_coding':
					id_ = attr['gene_id']
					name = attr['gene_name']
					if case_insensitive: name = name.upper()

					transcripts[id_] = name
					n_t += 1
			#if i >= 500000: break

	print "Parsed %d GFF/GTF lines, found %d valid gene entries, and %d valid transcript entries." %(i,n_g,n_t)

	# sanity check: do we have gene entries for all transcripts?
	all_transcripts = sorted(set(transcripts.keys()))
	for id_ in all_transcripts:
		assert id_ in genes
	print "There are %d genes with protein-coding transcripts." %(len(all_transcripts))

	known = 0
	for id_ in all_transcripts:
		if gene_status[id_] == 'KNOWN':
			known += 1
	print '%d / %d genes with protein-coding transcripts have gene_status "known".' %(known,len(all_transcripts))

	c = 0
	for id_ in all_transcripts:
		if id_.startswith('ENSGR'): c += 1
	print '%d / %d genes have IDs starting with "ENSGR" (genes in pseudoautosomal regions on X/Y).' %(c,len(all_transcripts))

	# associate names with chromosomes and IDs (not necessarily unique)
	name2chrom = {}
	name2id = {}
	for id_ in all_transcripts:
		name = transcripts[id_]
		chrom = chromosomes[id_]
		try:
			name2chrom[name].append(chrom)
			name2id[name].append(id_)
		except KeyError:
			name2chrom[name] = [chrom]
			name2id[name] = [id_]

	all_transcript_names = sorted(set(transcripts.values()))
	print "Among transcripts, found %d unique Ensembl IDs, and %d unique gene names" %(len(all_transcripts),len(all_transcript_names))

	mult = 0
	redundant_names = []
	redundant_ids = []
	for n,ids in name2id.iteritems():
		if len(ids) > 1:
			redundant_names.append(n)
			redundant_ids.append(ids)
			mult += (len(ids)-1)
	redundant = sorted(set(redundant_names))
	print mult,"cases of a gene name being associated with additional Ensembl IDs, affecting %d genes." %(len(redundant_names))
	#print redundant

	filtered_transcript_names = all_transcript_names[:]
	if exclude_chromosomes:
		e = 0
		filtered_transcript_names = []
		for name in all_transcript_names:
			chrom = sorted(set(name2chrom[name][:]))
			for c in exclude_chromosomes:
				try: chrom.remove(c)
				except ValueError: pass
			if chrom:
				filtered_transcript_names.append(name)
				name2chrom[name] = chrom
			else:
				e+=1
		if e > 0:
			print "Excluded %d genes based on excluded chromosomes!" %(e)
		

	with open(args.output_file,'w') as ofh:
		writer = csv.writer(ofh,dialect='excel-tab',lineterminator='\n',quoting=csv.QUOTE_NONE)
		for name in filtered_transcript_names:
			chroms = ','.join(sorted(set(name2chrom[name])))
			ids = ','.join(name2id[name])
			writer.writerow([name,chroms,ids])

	return 0

if __name__ == '__main__':
	return_code = main(read_args_from_cmdline())
	sys.exit(return_code)
