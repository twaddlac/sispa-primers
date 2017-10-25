#!/usr/bin/python

import sys
from itertools import combinations_with_replacement,product,islice,ifilter
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC as gc
from Bio.SeqUtils import lcc

for amp in combinations_with_replacement('ACTG', 4):
	
	amp = "".join(amp)
	print amp
	
	if any(nt*3 in amp for nt in 'ACTG'):
		continue
	# else:
	# 	print amp
		
