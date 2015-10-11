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

#import numpy as np
from gopca import common

def read_args_from_cmdline():

	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-b','--bootstrap-gopca-file',required=True)
	parser.add_argument('-i','--result-index',type=int,required=True)
	parser.add_argument('-o','--output_file',required=True)

	return parser.parse_args()

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()


	bootstrap_gopca_file = args.bootstrap_gopca_file
	result_index = args.result_index
	output_file = args.output_file

	assert os.path.isfile(bootstrap_gopca_file)
	assert result_index >= 0

	bootstrap_result = common.read_gopca_result(bootstrap_gopca_file)
	result = bootstrap_result.gopca_results[result_index]

	print 'Saving to file...', ; sys.stdout.flush()
	result.save(output_file)
	print 'done!'; sys.stdout.flush()
	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
