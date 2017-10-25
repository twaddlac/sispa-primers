#!/usr/bin/python

import sys
from itertools import product
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio.SeqUtils import lcc
from Bio.Seq import Seq
from primer3 import calcHairpin as hp
from primer3 import calcHomodimer as hd


# def generateStdAmps(p5,p7):
# 	amps = list()
# 	finalAmps = list()
# 	#print gc(p5[-11:]),gc(p7[-11:])
#
# 	for amp in product('ACTG', repeat=11):
#
# 		amp = "".join(amp)
# 		f = p5+amp
# 		f = f[-33:]
# 		r = p7+amp
# 		r = r[-33:]
#
# 		# no ^TT because of the T overhang on the adapter sequence
# 		if any(nt*3 in amp for nt in 'ACTG') or lcc.lcc_simp(amp) <= 1.8 or amp[:2] == "TT" or:
# 			continue
# 		# the generated amp needs to have 54% GC because the last 11bp of the adatpers have 45% GC and this will equate to 52% for seq primer
# 		elif 54 <= gc(amp) < 55 and 65 <= mt.Tm_NN(f) < 65:
# 			amps.append(amp)
#
# 	return amps

def findHairpin(sequence):
	longest = 0
	length = 1
	
	while length <= len(sequence)/2:
		comp = Seq.complement(Seq(sequence[:length]))
		if str(comp) in sequence[length:]:
			if length > longest:
				longest = length
		
		length += 1
	return longest
		
def findHomoDimer(sequence):
	length = 1
	longest = 0
	
	rev = sequence[::-1]
	while length <= len(sequence):
		if rev[:length] in sequence:
			if length > longest:
				longest = length
		length += 1
	
	return longest


def findBCSets():
	bcSets = dict()
	for bc in product('ACTG', repeat=7):
		barcode = "".join(bc)
		if lcc.lcc_simp(barcode) <= 1.8 or any(nt*3 in barcode for nt in 'ACTG'):
			continue
		elif len(bcSets) == 0:
			bcSets["set1"] = list()
			bcSets["set1"].append(barcode)
		else:
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
	return bcSets

def findPrimerSets(ampPrimers,p5,p7):
	count = 1
	seqPrimers = dict()
	for i in range(0,len(ampPrimers)):

		f = p5+ampPrimers[i]
		f = f[-33:]

		if 65 <= mt.Tm_NN(f) < 66:

			#print mt.Tm_NN(f),mt.Tm_NN(r),GC(f),GC(r)
			#print f,r,findHomoDimer(f),findHomoDimer(r)
			key = "set"+str(count)
			seqPrimers[key] = dict()
			seqPrimers[key]["seqAmp"] = f
			seqPrimers[key]["seqAmp_Tm"] = mt.Tm_NN(f)
			seqPrimers[key]["seqAmp_GC"] = GC(f)
			seqPrimers[key]["stdAmp"] = ampPrimers[i]
			count += 1

	return seqPrimers
	
	
def getValidSets(seqPrimers,bcSets,p5,p7):
	
	for seqSet,seqDict in seqPrimers.iteritems():
		for bcSet,bcs in bcSets.iteritems():
			bcCount = 0
			validSet = list()
			if len(bcs) >= 96:
				for bc in bcs:
					stdAmp = seqDict["stdAmp"] + bc
					
					if not any(nt*3 in stdAmp for nt in 'ACTG') and \
					62 <= mt.Tm_NN(stdAmp) < 64 and \
					GC(stdAmp) < 58 and \
					hp(stdAmp).tm < 50 and \
					hd(stdAmp).tm < 60:
						validSet.append(bc)
				#picking number 20 because 20, 22, 23 are all the same
				if len(validSet) >= 96:# and bcSet == "set20":
					print "primer_seq_seq barcode_set barcode stdamps stdamp_tm stdamp_gc stdamp_hairpin_tm stdamp_homodimer_dg	seq_primer seqprimer_tm seqprimer_gc seqprimer_hairpin_tm seqprimer_homodimer_dg"			
					for i in range(0,96):
						stdAmp = seqDict["stdAmp"]+validSet[i]
						seqAmp = seqDict["seqAmp"]
						print seqSet,bcSet,validSet[i],stdAmp,mt.Tm_NN(stdAmp),GC(stdAmp),hp(stdAmp).tm,str(hd(stdAmp).dg/1000),\
						seqAmp,mt.Tm_NN(seqAmp),GC(seqAmp),hp(seqAmp).tm,str(hd(seqAmp).dg/1000)
	

def main():
	
	#print findHairpin("AAACCCTTT")
	
	#print findHomoDimer("AAACCCTTT")

	# p5 is complement
	p5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	
	# p7 is reverse complement
	p7 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"

	# fuck jcvi primers
	# amps = list()
	# with open("jcvi-standard-amp-primers.txt",'rb') as f:
	# 	for line in f:
	# 		amps.append(line.rstrip())
	#
	# seqPrimers = findPrimerSets(amps,p5,p7)

	#seqPrimers = generateStdAmps(p5,p7)
	
	bcSets = findBCSets()
	
	getValidSets(seqPrimers,bcSets,p5,p7)

if __name__ == "__main__":
	main()
