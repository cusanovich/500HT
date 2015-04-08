#!/usr/bin/env python
import sys

if(len(sys.argv)!=3):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: tped2bed.py tped outbed\n\n')
	sys.exit(1)

tped = open(sys.argv[1],'r')
outbed = open(sys.argv[2],'w')
for line in tped:
	liner = line.strip().split()
	chrm = liner[0]
	end = liner[3]
	start = int(liner[3]) - 1
	rsid = liner[1]
	print >> outbed, "chr" + chrm + "\t" + str(start) + "\t" + end + "\t" + rsid + ("\t0\t.")

tped.close()
outbed.close()