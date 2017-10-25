#!/usr/bin/python

import sys
from itertools import product,islice,ifilter
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC as gc
from Bio.SeqUtils import lcc

amps = list()

p5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
p7 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"

#print gc(p5[-11:]),gc(p7[-11:])

for amp in product('ACTG', repeat=14):
	
	amp = "".join(amp)
	
	if any(nt*3 in amp for nt in 'ACTG') and lcc.lcc_simp(amp) >= 1.9:
		continue
	else:
		amps.append(amp)
		
for i in range(0,len(amps)):
	for x in range(0,len(amps)):
		if i == x:
			continue
		
		amp = amps[i]+amps[x]
		f = p5+amp
		f = f[-33:]
		r = p7+amp
		r = r[-33:]
				
		if not any(nt*3 in amp for nt in 'ACTG') and \
		i != x and \
		amps[i][::-1] != amps[x] and \
		amps[i] != amps[x][::-1] and \
		65 <= mt.Tm_NN(f) < 66 and \
		54 <= gc(amp) < 55:
			print f,r,mt.Tm_NN(f),mt.Tm_NN(r),gc(amp),amp

# for amp in product('ACTG',repeat=22):
# 	amp = "".join(amp)
# 	if any(nt*3 in amp for nt in 'ACTG'):
# 		continue
# 	else:
# 		print amp
#

