#!/usr/bin/env python2.7

import sys
import argparse
import cPickle as pickle

from scipy.io import savemat

from gopca import common

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-g','--gopca-file',required=True)
	parser.add_argument('-o','--output-file',required=True)
	parser.add_argument('--no-append-mat',action='store_true')

	return parser.parse_args()

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	gopca_file = args.gopca_file
	output_file = args.output_file
	append_mat = not args.no_append_mat

	result = None
	with open(gopca_file,'rb') as fh:
		result = pickle.load(fh)
	
	signatures = result.signatures
	labels = [sig.get_pretty_format() for sig in signatures]
	sig_genes = dict([sig.enrichment.term[0].replace(':','_'),sorted(sig.genes)] for sig in signatures)
	#print sig_genes
	samples = list(result.samples)
	S = result.S
	#common.write_expression(output_file,labels,samples,S)

	mat = {}
	mat['signature_labels'] = labels
	mat['samples'] = samples
	mat['S'] = S
	mat['signature_genes'] = sig_genes

	savemat(output_file,mat,appendmat=append_mat)

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
