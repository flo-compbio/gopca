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

import sys
import os
import argparse
import cPickle as pickle
import csv

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-g','--gopca-file',required=True)
	parser.add_argument('-o','--output-file',required=True)

	return parser.parse_args()

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	gopca_file = args.gopca_file
	output_file = args.output_file

	assert os.path.isfile(gopca_file)

	result = None
	with open(gopca_file,'rb') as fh:
		result = pickle.load(fh)
	
	signatures = result.signatures

	signatures = sorted(signatures,key=lambda sig:sig.term[3])

	with open(output_file,'w') as ofh:
		writer = csv.writer(ofh,dialect='excel-tab',lineterminator='\n',quoting=csv.QUOTE_NONE)
		# write header
		header = ['Term ID','Label','PC','No. of genes','Total no. of genes','Threshold','p-value','Max. fold enrichment','Genes']
		writer.writerow(header)
		for sig in signatures:
			gene_str = ','.join(sorted(sig.genes))
			data = [sig.term[0],sig.get_label(include_id=False),str(sig.pc),str(sig.k),str(sig.K),str(sig.n),'%.1e' %(sig.pval),'%.1f' %(sig.mfe),gene_str]
			writer.writerow(data)

	print 'Wrote %d signatures to "%s".' %(len(signatures),output_file)
	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
