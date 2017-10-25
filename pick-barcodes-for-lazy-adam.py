#!/usr/bin/python

import sys

bc = list()

with open(sys.argv[1],"rb") as f:
	
	for line in f:
		line = line.rsplit(" ")
		bc.append(line[1])
		
for i in range(0,len(bc)):
	valid = list()
	valid.append(bc[i])
	for x in range(i,len(bc)):
		if sum(ch1 != ch2 for ch1, ch2 in zip(bc[i],bc[x])) == 7 :
			valid.append(bc[x])
	if len(valid) < 8:
		break
	else:
		nt = dict()
		for x in range(0,7):
			nt[x] = dict()
			for nt in 'ACTG':
				nt[x][nt] = 0
		
		for x in range(0,7)
	