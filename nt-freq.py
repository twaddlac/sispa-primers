#!/usr/bin/python

import sys

nt = dict()
for i in range(0,7):
	nt[i] = dict()
	for x in 'ACTG':
		nt[i][x] = 0

with open(sys.argv[1],"rb") as f:
	for line in f:
		line.rstrip()
		for i in range(0,7):
			nt[i][line[i]] += 1

print "pos A C T G"
for k in nt:
	# print k,v['A'],v['C'],v['T'],v['G']
	print k,nt[k]['A'],nt[k]['C'],nt[k]['T'],nt[k]['G']