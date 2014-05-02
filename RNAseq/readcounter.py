#!/usr/bin/python
#This script was written by DC to count reads in multiplexed RNAseq samples 
#1/28/12, updated for flow cell ID on 2/4/12

import subprocess
import sys

if(len(sys.argv)!=2):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: readcounter.py FlowCell#\n\n')
	sys.exit(1)

fcno = sys.argv[1]
fcid = 'FlowCell' + sys.argv[1] + '/'
print fcid
fcidd = 'FlowCell' + sys.argv[1]
print fcidd
outter = '/mnt/lustre/data/users/cusanovich/500HTRNAseq/' + fcidd + 'b_500HT_RawReadCounts_byFindiv.txt'
print outter
outfile1 = open(outter,"w")

print 'Making count file...'
print >> outfile1, 'Findiv\tLane1\tIndex1\tReads1\tLane2\tIndex2\tReads2\tTotalReads'

##Generate a list of subjects with seq files
checker = 'ls /mnt/lustre/data/users/cusanovich/500HTRNAseq/' + fcid
q = subprocess.Popen(checker, shell=True,stdout=subprocess.PIPE)
findivs = q.communicate()[0]
finlist = findivs.split()
#print finlist[0:3]
i = 1

##This loop counts lines in each seq file to determine the number of reads per sample
for id in finlist:
	print 'Sample ' + str(i)
	currseqs = 'ls /mnt/lustre/data/users/cusanovich/500HTRNAseq/' + fcid + id + '/'
	seq = subprocess.Popen(currseqs, shell=True,stdout=subprocess.PIPE)
	seqs = seq.communicate()[0]
	seqlist = seqs.split()
	count = 0
	for sample in seqlist:
		target = 'zcat /mnt/lustre/data/users/cusanovich/500HTRNAseq/' + fcid + id + '/' + sample + ' | wc -l'
		print target
		p = subprocess.Popen(target, shell=True,stdout=subprocess.PIPE)
		pnum = p.communicate()[0].rstrip()
		samp = sample.split('.')
		if count == 0:
			lane1 = samp[0]
			index1 = samp[1]
			pnum1 = pnum
			count = 1
			continue
		if count == 1:
			lane2 = samp[0]
			index2 = samp[1]
			pnum2 = pnum
	print >> outfile1, id + "\t" + lane1 + "\t" + index1 + "\t" + str(int(pnum1)/4) + "\t" + lane2 + "\t" + index2 + "\t" + str(int(pnum2)/4) + "\t" + str(int(pnum1)/4 + int(pnum2)/4)
	print str(i) + ' of ' + str(len(finlist)) + ' be finished.'
	i += 1
outfile1.close()

infile = open(outter,"r")

line = infile.readline()
line = infile.readline()	#skip header
counts = {}
while line:
	fields = line.strip('\t').split()
	findiv = fields[0]
#	print findiv
	counts[ findiv ] = [ fields[1], fields[2], int(fields[3]), fields[4], fields[5], int(fields[6]), int(fields[7]) ]
	line = infile.readline()
totals = 0
lanes = [0,0,0,0,0,0,0,0]
indices = [0,0,0,0,0,0,0,0,0,0,0,0]
for key in counts.keys():
	totals = counts[key][6] + totals
	laneA = int(counts[key][0].split('_')[1]) - 1
	laneB = int(counts[key][3].split('_')[1]) - 1
	indexA = int(counts[key][1].split('_')[1]) - 1
	indexB = int(counts[key][4].split('_')[1]) - 1
	lanes[laneA] = lanes[laneA] + counts[key][2]
	indices[indexA] = indices[indexA] + counts[key][2]
	lanes[laneB] = lanes[laneA] + counts[key][5]
	indices[indexB] = indices[indexA] + counts[key][5]
infile.close()
outfile = open(outter,"a")
print >> outfile, "Total Reads = " + str(totals)
print >> outfile, "\n"
print >> outfile, "Lane Totals"
print >> outfile, "Lane1\tLane2\tLane3\tLane4\tLane5\tLane6\tLane7\tLane8"
print >> outfile, "\t".join(map(str,lanes))
print >> outfile, "\n"
print >> outfile, "Index Totals"
print >> outfile, "Index1\tIndex2\tIndex3\tIndex4\tIndex5\tIndex6\tIndex7\tIndex8\tIndex9\tIndex10\tIndex11\tIndex12"
print >> outfile, "\t".join(map(str,indices))
outfile.close()
