#!/usr/bin/python

import sys
from itertools import product

# BC004CG
sa1 = "CGTAGTACACTCTAGAGCACTA"

# BC009CG
#sa2 = "CGAGCTCTATACGTGTAGTCTC"

bcSet = []
names = []

i = 1

for bc in product('ACTG', repeat=6):
	match = True
	barcode = sa1+"".join(bc)
	if len(bcSet) == 0:
		bcSet.append(barcode)
		names.append(str(i)+"".join(bc))
	else:
		for x in range(0,len(bcSet)):
			print bcSet[x],barcode
			if sum(ch1 != ch2 for ch1, ch2 in zip(bcSet[x],barcode)) < 3:
				match = False
				break
				#print sum(ch1 != ch2 for ch1, ch2 in zip(bcSet[x],barcode))
		
		if match == True:
			bcSet.append(barcode)
			names.append(str(i)+"".join(bc))
			
	print len(bcSet)
	i += 1
	

with open ("barcodes.fa","w") as f:
	for i in range(0,len(bcSet)):
		f.write(">"+str(i)+"\n"+bcSet[i]+"\n")

# ham = dict()
#
# for i in range(len(a1)):
# 	for x in range (i+1,len(a1)):
# 		ham[i][x] = 0

# print "\t"+"\t".join(names)
# for i in range(len(a1)):
# 	print names[i],
# 	for x in range (i+1,len(a1)):
# 		print str(sum(ch1 != ch2 for ch1, ch2 in zip(a1[i], a1[x])))+"\t",
# 	print
		#ham[i][x] = sum(ch1 != ch2 for ch1, ch2 in zip(a1[i], a1[x]))
		
