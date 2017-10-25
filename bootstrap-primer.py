#!/usr/bin/python

# given a sequence, it shuffles the bases and tests it to see if its better than the previous

import sys
from itertools import *
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC as gc
from Bio.SeqUtils import lcc
from Bio.Seq import Seq
from primer3 import calcHairpin as hp
from primer3 import calcHomodimer as hd
import time
import random
import numpy as np


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
	
	finalBCSets = dict()
	for bc,bcList in bcSets.iteritems():
		if len(bcList) >= 96:
			finalBCSets[bc] = bcList
	
	return finalBCSets


def testBarcodes(bcSets,primer):
	
	amp = primer[-22:]

	validSets = dict()
	
	for bcSet,bcList in bcSets.iteritems():
		worksList = list()
		for bc in bcList:
			seq = amp+bc
			
			if 60 < mt.Tm_NN(seq) < 65 and hp(seq).tm <= 42 and hd(seq).dg/1000 > 2:
				worksList.append(bc)
				#print hd(seq).dg/1000
				# print  mt.Tm_NN(seq),hp(seq).tm,hd(seq).dg/1000
		if len(worksList) >= 96:
			validSets[bcSet] = worksList
	
	return validSets

def main():
	
	#orig seq. had to replace last G with T to get to 52 pct
	# this was taken from JCVI's list, but as long as it has 52 pct GC it doesn't matter
	# seqPrimer = "TCTTCCGATCTCGTCGCTACGCGTCATCAGTAG"
	
	seqPrimer = "TCTTCCGATCTCGTCGCTACGCGTCATCAGTAT"
	
	lastTm = mt.Tm_NN(seqPrimer)
	lastHP = hp(seqPrimer).tm
	lastHD = hd(seqPrimer).dg/1000
	lastgc = gc(seqPrimer)
	lastSeqPrimer = seqPrimer
	validSets = list()
	
	bcSets = findBCSets()
	
	count = 0
	
	while count <= 1000000:
		# shuffles the sequences while retaining the GC%
		amp = np.random.permutation(list(seqPrimer[-22:]))
		amp = "".join(amp)
		if any(nt*3 in amp for nt in 'ACTG'):
			continue
	
		amp = "".join(amp)
		amp = seqPrimer[:11]+amp
		newTm = mt.Tm_NN(amp)
		newHP = hp(amp).tm
		newHD = hd(amp).dg/1000
		
		if newTm < 65:
			continue
		elif lastHP >= newHP and lastHD < newHD:
			
			print amp,lastTm,newTm,lastHP,newHP,lastHD,newHD
			
			validSets = testBarcodes(bcSets,amp)
			
			if len(validSets) > 0:
						
				lastHP = newHP
				lastHD = newHD
				lastTm = newTm
				lastSeqPrimer = amp

		count += 1

	print amp,lastTm,newTm,lastHP,newHP,lastHD,newHD
		
		
	
if __name__ == "__main__":
	main()