#!/usr/bin/python
## This script was written by DC to run fastqc on sequencing files
## 1/27/12

import sys, os, subprocess

if(len(sys.argv)!=2):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: fastqc.py FlowCell#\n\n')
	sys.exit(1)

fcid = sys.argv[1]
print fcid

outfile = open("/mnt/lustre/home/caliskan/fc" + fcid + "fastqc.sh","w")


print >> outfile, "#!/bin/bash"


## Generate a list of subjects with seq files
q = subprocess.Popen('ls /mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell' + fcid + '/', shell=True,stdout=subprocess.PIPE)
findivs = q.communicate()[0]
finlist = findivs.split()
## This loop identifies all of the samples sequenced to generate a bash script for running fastqc
for id in finlist:
	currseqs = 'ls /mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell' + fcid + '/' + id + '/*.sequence.txt.gz'
#	print currseqs
	seq = subprocess.Popen(currseqs, shell=True,stdout=subprocess.PIPE)
	seqs = seq.communicate()[0]
	seqlist = seqs.split()
	count = 0
#	print seqlist[1]
	for sample in seqlist:
		print >> outfile, "/mnt/lustre/home/jroux/bin/FastQC/fastqc --outdir=/mnt/lustre/home/caliskan/FASTQC/FlowCell" + fcid + "/ " + sample
outfile.close()