#!/usr/bin/python
## This script was written by DC to reorganize multiplexed RNAseq samples - it uses the "SampleDirectories.csv" file to identify which subdirectories different
## samples ended up in and then copies the sequence.txt.gz files to a directory of your choice.
## 1/27/12

import sys, os, subprocess

if(len(sys.argv)!=2):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: fastqc.py FlowCell#\n\n')
	sys.exit(1)

fcid = sys.argv[1]
print fcid

outfile = open("/mnt/lustre/data/users/cusanovich/RNAseq_scripts/fc" + fcid + "fastqc.sh","w")


print >> outfile, "#!/bin/bash"


## Generate a list of subjects with seq files
q = subprocess.Popen('ls /mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell' + fcid + '/', shell=True,stdout=subprocess.PIPE)
findivs = q.communicate()[0]
finlist = findivs.split()
#print finlist[0:3]
#i = 1
## This loop identifies all of the samples sequenced to generate a bash script for running fastqc
for id in finlist:
	if id  == "Saved":
		continue
	currseqs = 'ls /mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell' + fcid + '/' + id + '/*.sequence.txt.gz'
#	print currseqs
	seq = subprocess.Popen(currseqs, shell=True,stdout=subprocess.PIPE)
	seqs = seq.communicate()[0]
	seqlist = seqs.split()
	count = 0
#	print seqlist[1]
	for sample in seqlist:
#		print "fastqc --outdir=/mnt/lustre/data/users/cusanovich/500HTRNAseq/FASTQC/FlowCell" + fcid + "/ " + sample
		print >> outfile, "fastqc --outdir=/mnt/lustre/data/users/cusanovich/500HTRNAseq/FASTQC/FlowCell" + fcid + "/ " + sample
outfile.close()
