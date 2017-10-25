#!/usr/bin/python

import sys
from itertools import product
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio.SeqUtils import lcc
from Bio.Seq import Seq
# read 1 adapter
# always read in the forward direction
p5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTA"

# read 2 adapter
# this needs to be reverse complemented as well as the sequence 5' of it.
p7 = Seq.reverse_complement(Seq("ATCTCGTATGCCGTCTTCTGCTTG"))

amps = list()
with open("jcvi-standard-amp-primers.txt",'rb') as f:
	for line in f:
		amps.append(line.rstrip())

rcAmps = list()
for i in range(0,len(amps)):
	rcAmps.append(Seq.reverse_complement(Seq(amps[i])))



# for amp in amps:
# 	ampPrimer = adapt+amp
# 	print ampPrimer[-32:],mt.Tm_NN(ampPrimer[-32:]),GC(ampPrimer[-32:])

bcSets = dict()

for bc in product('ACTG', repeat=7):
	barcode = "".join(bc)
	if lcc.lcc_simp(barcode) <= 1.8 or any(nt*3 in barcode for nt in 'ACTG'):
		continue
	elif len(bcSets) == 0:
		#print barcode,lcc.lcc_simp(barcode)
		bcSets["set1"] = list()
		bcSets["set1"].append(barcode)
	else:
		#print barcode,lcc.lcc_simp(barcode)
		setsTried = 0
		for k,v in bcSets.iteritems():
			match = True
			setsTried += 1
			for i in range(0,len(v)):
				if sum(ch1 != ch2 for ch1, ch2 in zip(v[i],barcode)) < 3:
					match = False
					break

			if match == True:
				v.append(barcode)
				setsTried = 0
				break
			else:
				if setsTried == len(bcSets):
					bcSets["set"+str(setsTried+1)] = list()
					bcSets["set"+str(setsTried+1)].append(barcode)
					break
				elif setsTried == len(bcSets):
					bcSets["set"+str(setsTried)] = list()
					bcSets["set"+str(setsTried)].append(barcode)
					break


print "set\tamp primer\tbarcode\tbarcode complexity\tprimer Tm\tprimer GC\tampPrimeruencing primer\tampPrimeruencing primer Tm\tampPrimeruencing primer GC"

# finalDict = dict()

for k,v in bcSets.iteritems():
	for amp in amps:
		# finalDict[k+amp] = list()
		for i in range(0,len(v)):
			ampPrimer = amp+v[i]
			seqPrimer = adapt+amp
			if mt.Tm_NN(seqPrimer[-33:]) >= 65 and 49 <= GC(seqPrimer[-33:]) <= 55 and mt.Tm_NN(ampPrimer) <= 64 and not any(nt*3 in ampPrimer for nt in 'ACTG'):
				print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(k,v[i],lcc.lcc_simp(v[i]),ampPrimer,mt.Tm_NN(ampPrimer),
				GC(ampPrimer),seqPrimer[-33:],mt.Tm_NN(seqPrimer[-33:]),GC(seqPrimer[-33:]))
				# finalDict[k+amp].append(v[i]+"\t"+str(lcc.lcc_simp(v[i]))+"\t"+ampPrimer+"\t"+str(mt.Tm_NN(ampPrimer))+"\t"+
				# str(GC(ampPrimer))+"\t"+seqPrimer[-33:]+"\t"+str(mt.Tm_NN(seqPrimer[-33:]))+"\t"+str(GC(seqPrimer[-33:])))

# seqPrimerTm = 100
# tmList = list()
# for k,v in finalDict.iteritems():
# 	if len(v) >= 96:
# 		line = v[0].split("\t")
# 		if abs(float(line[6]) - 65) < seqPrimerTm:
# 			seqPrimerTm = abs(float(line[6]) - 65)
# 			tmList = list()
#
# 			for i in range(0,len(v)):
# 				if len(tmList) <= 96:
# 					tmList.append(v[i])
# 					continue
# 				else:
# 					for x in range(0,len(tmList)):
# 						tm = tmList[x].split("\t")
# 						if float(line[6]) <= float(tm[6]):
# 							tmList.insert(x,v[i])
# 							last = tmList.pop
#
# for x in tmList:
# 	print x
							
							
							
							
				
				
				
				
				