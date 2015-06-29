#!/usr/bin/env python

import sys
import os
import argparse
import csv
import gzip

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-f','--gene2acc-file',required=True)
	parser.add_argument('-o','--output-file',required=True)

	return parser.parse_args()

def open_plain_or_gzip(fn):
	try:
		gzip.open(fn).next()
		return gzip.open(fn)
	except IOError:
		return open(fn)

def read_gene2acc(fn):
	gene2acc = {}
	with open_plain_or_gzip(fn) as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		reader.next() # skip header
		for i,l in enumerate(reader):
			if (i % 100000) == 0:
				print i, ; sys.stdout.flush()

			id_ = int(l[1])
			symbol = l[15]

			try:
				gene2acc[id_].append(symbol)
			except:
				gene2acc[id_] = [symbol]
			#print l[0],l[15]
	print

	# make sure all EntrezIDs map to a unique gene symbol
	n = len(gene2acc.keys())
	for k,v in gene2acc.iteritems():
		symbols = sorted(set(v))
		assert len(symbols) == 1
		gene2acc[k] = symbols[0]

	all_symbols = sorted(set(gene2acc.values()))
	m = len(all_symbols)

	print "Found %d Entrez Gene IDs associated with %d gene symbols." %(n,m)
	sys.stdout.flush()
	return gene2acc

def write_entrez2gene(ofn,entrez2gene):
	with open(ofn,'w') as ofh:
		writer = csv.writer(ofh,dialect='excel-tab',lineterminator='\n')
		for k,v in entrez2gene.iteritems():
			writer.writerow([k,v])

def main():

	args = read_args_from_cmdline()
	entrez2gene = read_gene2acc(args.gene2acc_file)
	write_entrez2gene(args.output_file,entrez2gene)
	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
