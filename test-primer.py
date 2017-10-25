#!/usr/bin/python

# tested the primers based on what boostrap-primer.py output. Can be put into one script, just didn't 

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

# this was taken from 'boostrap-primer.py' 
seqPrimer = "TCTTCCGATCTGCCGCCTTACACTGACTTGAGT"

bcSets = findBCSets()

validSets = testBarcodes(bcSets,seqPrimer)

print "bcSet bc seqPrimer seqPrimerGC seqPrimerTm seqPrimerHP seqPrimerHD stdAmp stdAmpGC stdAmpTm stdAmpHP stdAmpHD"
for k,v in validSets.iteritems():
	for bc in v:
		print k,bc,seqPrimer,gc(seqPrimer),mt.Tm_NN(seqPrimer),hp(seqPrimer).tm,hd(seqPrimer).dg/10000,seqPrimer[-22:]+bc,gc(seqPrimer[-22:]+bc),mt.Tm_NN(seqPrimer[-22:]+bc),hp(seqPrimer[-22:]+bc).tm,hd(seqPrimer[-22:]+bc).dg/1000
		