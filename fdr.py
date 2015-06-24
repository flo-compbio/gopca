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

import numpy as np

def fdr_bh(p_values,q):
	# implements Benjamini Hochberg procedure
	a = np.argsort(p_values)
	n = a.size
	k = 0

	for i in range(n):
		if p_values[a[i]] <= q*((i+1)/float(n)):
			k = i+1

	crit_p = 0.0
	if k > 0:
		crit_p = p_values[a[k-1]]

	return k,crit_p


def fdr_bh_general(p_values,q):
	# implements Benjamini Hochberg procedure that guarantees
	# FDR control under arbitrary test dependencies
	n = p_values.size

	q_adj = 0
	for i in range(n):
		q_adj += (1/float(i+1))
	q_adj = q/q_adj

	return fdr_bh(p_values,q_adj)
