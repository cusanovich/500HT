#!/usr/bin/python
#This fragment was written by DC to count reads in multiplexed RNAseq samples by lane and by index
#1/29/12

import subprocess,sys

if(len(sys.argv)!=2):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: readcounter.py FlowCell#\n\n')
	sys.exit(1)

fcno = sys.argv[1]
fcidd = 'FlowCell' + sys.argv[1]
print fcidd
counter = '/mnt/lustre/data/users/cusanovich/500HTRNAseq/' + fcidd + '_500HT_RawReadCounts_byFindiv.txt'

infile = open(counter,"r")

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
outfile = open(counter,"a")
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
#totals = counts.values()[6]
#print totals
#laneA = []
#for sample in counts:
#	print sample
#	info = counts.values()
#	for datum in info:
#		laneA.append(datum[0])
#	print temp
#	laneA = temp[1]
