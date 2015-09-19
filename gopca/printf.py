from __future__ import print_function

import sys

def printf(s, end=' \n', fh=sys.stdout):
	print(s,end=end,file=fh)
