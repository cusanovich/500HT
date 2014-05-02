#!/usr/bin/python
# Script written by DC to collect usable reads from those with bad
# indices. It was inspired by the perl script written by Russ Bainer
# for the same purpose. It differs in a few respects:
# (a) it only allows 1 mismatch between the sequenced index and the
# intended index;
# (b) it saves the usable reads in a file called *.saved.sequence.txt.gz
# rather than adding them to the original seq file.
# 2/08/12, edited 2/09/12 to run faster
# /mnt/lustre/data/share/SolexaSequencer/BaseCalls/120907_SN_0795_0145_BC169DACXX/Demultiplexed/unknown/GERALD_13-09-2012_root/

import sys
import os
import subprocess
import gzip


if(len(sys.argv) != 3):
    sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\n \
    Usage: gleaner.py UnknownLocation FlowCell#\n\n')
    sys.exit(1)

indir = sys.argv[1]
fcno = sys.argv[2]
fcid = 'FlowCell' + sys.argv[2] + '/'
#print fcid
outdir = '/mnt/lustre/data/users/cusanovich/500HTRNAseq/' + 'Saved/' + fcid

indexes = open("/mnt/lustre/data/users/cusanovich/RNAseq_scripts/indexes.txt", "r")
indices = {}
for line in indexes:
    currindex = line.strip().split()
    indices[currindex[0]] = currindex[1]

checker = 'ls ' + indir + '*sequence.txt.gz'
q = subprocess.Popen(checker, shell=True, stdout=subprocess.PIPE)
lanes = q.communicate()[0]
lanelist = lanes.split()

for lane in lanelist:
#    if lane.split('/')[-1] != 's_7_sequence.txt.gz' and  lane.split('/')[-1] != 's_8_sequence.txt.gz':
#        continue
    currlane = gzip.open(lane, "r")
    line = currlane.readline()
    laned = lane.strip().split("/")[-1].split("_")[1]
    indexedfiles = {}
    for item in indices:
        currout = outdir + 'lane_' + laned + '.index_' + indices[item] + '.saved.sequence.txt.gz'
        indexedfiles[item] = gzip.open(currout, "w")
    while line:
        if line[0:1] == "@":
            liner = line.strip().split('#')[1]
            tag = liner.split('/')[0]
            missing = 6
            winner = 'none'
            for index in indices:
                currindex = index
                mismatch = 0
                for base in range(0, 6):
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
