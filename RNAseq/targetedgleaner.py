#!/usr/bin/python
# This script was written by DC to collect usable reads from those with bad indices. It was inspired by the perl script written by Russ Bainer for the same purpose.
# It differs in a few respects: (a) it only allows 1 mismatch between the sequenced index and the intended index; (b) it saves the usable reads in a file called 
# *.saved.sequence.txt.gz rather than adding them to the original seq file.
# 2/08/12, edited 2/09/12 to run faster
# FC1: /mnt/lustre/data/share/SolexaSequencer/BaseCalls/120113_SN_0795_0094_AC08WMACXX/Demultiplexed/unknown/GERALD_20-01-2012_root/
# FC2: /mnt/lustre/data/share/SolexaSequencer/BaseCalls/120123_SN_0795_0096_AD0ELLACXX/Demultiplexed/unknown/GERALD_31-01-2012_root/
# FC3: /mnt/lustre/data/share/SolexaSequencer/BaseCalls/120308_SN_0795_0103_AD0KKBACXX/Demultiplexed/unknown/GERALD_20-03-2012_root/

import sys, os, subprocess, gzip


if(len(sys.argv)!=5):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: targetedgleaner.py UnknownLocation FlowCell# TargetIndex TargetLane\n\n')
	sys.exit(1)

indir = sys.argv[1]
fcno = sys.argv[2]
fcid = 'FlowCell' + sys.argv[2] + '/'
#print fcid
outdir = '/mnt/lustre/data/users/cusanovich/500HTRNAseq/Saved/' + fcid
targetindex = sys.argv[3]
targetlane = "s_" + sys.argv[4] + "_sequence.txt.gz"
print "Target index is " + str(targetindex)
print "Target lane is " + targetlane

## Grabs index sequences from a file in the scripts directory
indexes = open("/mnt/lustre/data/users/cusanovich/RNAseq_scripts/indexes.txt","r")
indices={}
for line in indexes:
	if line.strip().split("\t")[1] == targetindex:
		currindex = line.strip().split()
		indices[currindex[0]] = currindex[1]
print indices

#checker = 'ls ' + indir + '*sequence.txt.gz'
#q = subprocess.Popen(checker, shell=True,stdout=subprocess.PIPE)
#lanes = q.communicate()[0]
#lanelist = lanes.split()

## Kept the for loop for consistency between gleaner.py and targetedgleaner.py, but basically just opens seq file and looks for tags matching (w/in 1 mismatch) the tag of interest
for count in 'a':
	lane = indir + targetlane
	print "Lane = " + str(lane)
	currlane = gzip.open(lane, "r")
	line = currlane.readline()
	laned = lane.strip().split("/")[-1].split("_")[1]
	indexedfiles = {}
	for item in indices:
		currout = outdir + 'lane_' + laned + '.index_' + indices[item] + '.saved.sequence.txt.gz'
		indexedfiles[item] = gzip.open(currout,"w")
		print indexedfiles
	while line:
		if line[0:1] == "@":
			liner = line.strip().split('#')[1]
			tag = liner.split('/')[0]
			missing = 6
			winner = 'none'
			for index in indices:
				currindex = index
				mismatch = 0
				for base in range(0,6):
					if tag[base] != index[base]:
						mismatch += 1
				if mismatch < missing:
					winner = index
					missing = mismatch
					continue
				if mismatch == missing:
					winner = 'none'
			if missing < 2 and winner != 'none':
				outter = indexedfiles[winner]
				print >> outter, line.strip()
				print >> outter, currlane.readline().strip()
				print >> outter, currlane.readline().strip()
				print >> outter, currlane.readline().strip()
		line = currlane.readline()
	for item in indexedfiles:
		indexedfiles[item].close()
#zipper = 'gzip -r ' + outdir
#subprocess.Popen(zipper, shell=True)