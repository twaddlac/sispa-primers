#!/usr/bin/python

import sys
from itertools import product,combinations_with_replacement
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC as gc
from Bio.SeqUtils import lcc
from Bio.Seq import Seq
from primer3 import calcHairpin as hp
from primer3 import calcHomodimer as hd
import time

p5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
p7 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"

start = time.clock()

count = 0
#with open("stdamps.txt","w") as out:
	
for amp in product('ATGC',repeat=22):
	
	count += 1

	amp = "".join(amp)
	print amp
	if amp[:3] == "AAA":
		continue
	
	print amp
	# amp = p5+amp
	# amp = amp[-33:]

	# if not any(nt*3 in amp for nt in 'ACTG'):

	if count%10000000 == 0:
		print count,str(time.clock()-start)
	
print count	
		# if any(nt*3 in amp for nt in 'ACTG'):
		# 	continue
		#
		# if lcc.lcc_simp(amp) >= 1.9 and 54 <= gc(amp) < 55 and 65 <= mt.Tm_NN(amp) < 66 \
		# and hp(amp).tm < 50 and hd(amp).dg/1000 > 0:
		# 	out.write(amp)
		# else:
		# 	continue
