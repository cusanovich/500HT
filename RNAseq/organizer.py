#!/usr/bin/python
#This script was written by DC to reorganize multiplexed RNAseq samples - it uses the "SampleDirectories.csv" file to identify which subdirectories different samples ended up in
# and then copies the sequence.txt.gz files to a directory of your choice.
#1/27/12
#/mnt/lustre/data/share/SolexaSequencer/BaseCalls/120113_SN_0795_0094_AC08WMACXX/

import sys
import os
import csv
import subprocess
#from subprocess import call


if(len(sys.argv)!=3):
    sys.stderr.write('Error!! Only supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: organizer.py InDirectory FlowCell#\n\n')
    sys.exit(1)

def canibeint(findiv):
    """Return true if value can be coerced to int.

    Source:
    http://stackoverflow.com/questions/1265665/python-check-if-a-string-represents-an-int-without-using-try-except

    """
    try:
        int(findiv)
        return True
    except ValueError:
        return False


indir = sys.argv[1]
outdir = "/mnt/lustre/data/users/cusanovich/500HTRNAseq/"
fcno = sys.argv[2]

## Gerald inserts a directory level between the index and the sequences.  This just figures out the name of that directory.
interdircheck = "ls " + sys.argv[1] + "Demultiplexed/001/"
interdircall = subprocess.Popen(interdircheck, shell=True,stdout=subprocess.PIPE)
interdircom = interdircall.communicate()[0]
interdir = "/" + str(interdircom.split()[0])


indexes = open("/mnt/lustre/data/users/cusanovich/RNAseq_scripts/indexes.txt","r")
outfile = open("/mnt/lustre/data/users/cusanovich/RNAseq_scripts/fc" + fcno + "organizer.sh","w")

print >> outfile, "#!/bin/bash"

## Create a dictionary of indices
indices={}
for line in indexes:
    currindex = line.strip().split()
    indices[currindex[0]] = currindex[1]

## Opens sample record (file called 'SamplesDirectories.csv') and make a bash script to copy all of the sequence files based on the sample directory information
samples = file(indir + 'Demultiplexed/SamplesDirectories.csv')
samples.next()
count = 1 
for row in csv.reader(samples):
    if not canibeint(row[2]):
        continue
    lane = row[1]
    findiv = row[2]
    index = row[4]
    indno = indices[index]
    subdir = row[9]
    newdir = "mkdir " + outdir + "FlowCell" + fcno + "/" + findiv
    subprocess.call(newdir, shell=True)
    cpfile = "cp " + indir + "Demultiplexed/" + subdir + interdir + "/s_" + lane + "_sequence.txt.gz " + outdir + "FlowCell" + fcno + "/" + findiv + "/lane_" + lane + ".index_" + indno + ".sequence.txt.gz"
    counter = "echo Starting No. " + str(count) + "..."
    print >> outfile, counter
    count += 1
    print >> outfile, cpfile
print >> outfile, "echo \"Yeah, boy!\""
outfile.close()
copying = "bash /mnt/lustre/data/users/cusanovich/RNAseq_scripts/fc" + fcno + "organizer.sh"
subprocess.call(copying, shell=True)