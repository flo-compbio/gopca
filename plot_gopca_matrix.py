import sys
import argparse
import cPickle as pickle

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc

"""
if __name__ == '__main__' and __package__ is None:
	# allow explicit relative imports in executable script
	# source: http://stackoverflow.com/a/6655098
	parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	sys.path.insert(1, parent_dir)
	import gopca
	__package__ = 'gopca'
"""

#import gopca
from gopca import common
from gopca.tools import misc
from gopca.go_pca_objects import mHGTermResultWithPC

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-e','--expression-file',required=True)
	parser.add_argument('-s','--signature-file',required=True)
	parser.add_argument('-g','--go-pickle-file',required=True)

	return parser.parse_args()

def main(args):

	expression_file = args.expression_file
	signature_file = args.signature_file
	go_pickle_file = args.go_pickle_file

	# read expression data
	genes,samples,E = common.read_expression(expression_file)

	# read GO pickle
	GO = None
	with open(go_pickle_file) as fh:
		GO = pickle.load(fh)

	# read signatures
	signatures = None
	with open(signature_file) as fh:
		signatures = pickle.load(fh)

	# constructing signatures
	S = np.row_stack([common.get_signature_expression(genes,E,sig.genes) for sig in signatures])
	print S.shape; sys.stdout.flush()
	labels = [common.get_signature_label(GO,sig) for sig in signatures]
	q,n = S.shape
	#print m

	# clustering of rows
	distxy = squareform(pdist(S, metric='correlation'))
	R = dendrogram(linkage(distxy, method='average'),no_plot=True)
	order_rows = np.int64([int(l) for l in R['ivl']])
	S = S[order_rows,:]
	labels = [labels[idx] for idx in order_rows]

	# clustering of columns
	distxy = squareform(pdist(S.T, metric='euclidean'))
	R = dendrogram(linkage(distxy, method='average'),no_plot=True)
	order_cols = np.int64([int(l) for l in R['ivl']])
	S = S[:,order_cols]

	# plotting
	rc('font',family='serif',size=32)
	rc('figure',figsize=(18,16))

	# plotting
	plt.imshow(S,interpolation='none',aspect='auto',vmin=-3,vmax=3)
	plt.yticks(np.arange(q),labels,size=14)
	plt.xlabel('Samples (n=%d)' %(n),size='small')
	plt.show()

	rc('font',family='serif',size=32)
	rc('figure',figsize=(20,16))

	plt.imshow(np.corrcoef(S),interpolation='none',aspect='auto',vmin=-1,vmax=1)
	cbar = plt.colorbar(shrink=0.5)
	cbar.ax.tick_params(labelsize=20)
	cbar.set_label('Pearson correlation',size=20)
	plt.yticks(np.arange(q),labels,size=14)
	plt.xticks(np.arange(q),())

	plt.show()

if __name__ == '__main__':
	return_code = main(read_args_from_cmdline())
	sys.exit(return_code)
