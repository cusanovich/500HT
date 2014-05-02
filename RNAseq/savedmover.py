#!/usr/bin/python
#This script was written by DC to write a bash script for reorganizing the saved sequence files 
#2/10/12

import subprocess,sys

if(len(sys.argv)!=2):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: savedmover.py FlowCell#\n\n')
	sys.exit(1)

fcno = sys.argv[1]
fcid = 'FlowCell' + sys.argv[1] + '/'
print fcid
fcidd = 'FlowCell' + sys.argv[1]
print fcidd
sloc = '/mnt/lustre/data/users/cusanovich/500HTRNAseq/Saved/' + fcid
rloc = '/mnt/lustre/data/users/cusanovich/500HTRNAseq/' + fcid
guide = '/mnt/lustre/data/users/cusanovich/500HTRNAseq/' + fcidd + '_500HT_RawReadCounts_byFindiv.txt'
guided = open(guide,"r")
outter = open("fc" + fcno + "savedmover.sh","w")
print >> outter, "#!/bin/bash"

line = guided.readline()
line = guided.readline()
while line:
	sample = line.strip().split()
	if sample[0] == 'Total':
		break
	findiv = sample[0]
	laneA = sample[1]
	indexA = sample[2]
	laneB = sample[4]
	indexB = sample[5]
	moverA = "cp " + sloc + laneA + "." + indexA + ".saved.sequence.txt.gz " + rloc + findiv + "/"
	moverB = "cp " + sloc + laneB + "." + indexB + ".saved.sequence.txt.gz " + rloc + findiv + "/"
	print >> outter, moverA
	print >> outter, moverB
	line = guided.readline()
outter.close()
copying = "bash /mnt/lustre/data/users/cusanovich/RNAseq_scripts/fc" + fcno + "savedmover.sh"
subprocess.call(copying, shell=True)