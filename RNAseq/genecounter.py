#!/usr/bin/python
# This script was written by DC to generate bed files for all the mapped reads
# 4/10/12

import os
import subprocess

# Step 1: generate findiv list
# Step 2: generate lane id list
# Step 3: merge quality.bam and saved.quality.bam files and generate bed
# Step 4: merge junction.quality.bam and saved.junction.quality.bam files, convert coordinates and generate bed
# Step 5: join bed files and re-sort

def findiv_grabber(fcid):
    """Return the names of the subdirectories (findivs) in a given directory (flow cell)."""
    checker = 'ls /mnt/lustre/data/users/cusanovich/500HTRNAseq/' + fcid
    q = subprocess.Popen(checker, shell=True,stdout=subprocess.PIPE)
    findivs = q.communicate()[0]
    finlist = findivs.split()
    return finlist

def lane_grabber(flowcell,findiv):
    """Return the the index and lane id for each replicate (that has a ".counts" file) in a sample directory.

    Return ids in the format: "path/of/interest/lane_X.index_Y"

    """
    currseqs = 'ls /mnt/lustre/data/users/cusanovich/500HTRNAseq/' + flowcell + findiv + '/*.counts.txt'
    counter = subprocess.Popen(currseqs, shell=True,stdout=subprocess.PIPE)
    counters = counter.communicate()[0]
    countslist = counters.split()
    countslist = [counted.rstrip('.counts.txt') for counted in countslist]
    return countslist

if(__name__=='__main__'):
    x = 0
    outdir = '/mnt/lustre/data/users/cusanovich/RNAseq_scripts/bedmakerOE/'
    for j in range(7,8):
        fcno = str(j)
        fcid = 'FlowCell' + fcno + '/'
        print fcid
        finlist = findiv_grabber(fcid)
        print finlist[0]

        samples = []
        for id in finlist:#['163802','168851']:#
            countslist = lane_grabber(fcid,id)
            samples.extend(countslist)
        print samples[0:11]

        for sample in samples:
            if os.path.isfile(sample + '.genecounts.txt'):
                print 'Skipped ' + "/".join(sample.split('/')[-2:])
                continue
#            if x >= 32:
#                print 'Broke for quota.'
#                break
#                      sort -k1,1 -k2,2n -k3,3 ' + sample + '.bed > ' + sample + '.sorted.bed; \
#                      rm ' + sample + '.bed; \
            merger = 'echo " \
                      /mnt/lustre/home/cusanovich/Programs/samtools/samtools merge ' + sample + '.quality.merged.bam ' + sample + '.quality.sort.bam ' + sample + '.saved.quality.sort.bam; \
                      /mnt/lustre/home/cusanovich/Programs/BEDTools/bin/bamToBed -i ' + sample + '.quality.merged.bam > ' + sample + '.bed; \
                      /mnt/lustre/home/cusanovich/Programs/samtools/samtools merge ' + sample + '.junction.quality.merged.bam ' + sample + '.junction.quality.sort.bam ' + sample + '.saved.junction.quality.sort.bam; \
                      /mnt/lustre/home/cusanovich/Programs/samtools/samtools view ' + sample + '.junction.quality.merged.bam > ' + sample + '.junction.quality.merged.sam; \
                      python /mnt/lustre/data/users/cusanovich/RNAseq_scripts/junctionreformatter.py ' + sample + '.junction.quality.merged.sam ' + sample + '.junction.bed; \
                      cat ' + sample + '.junction.bed >> ' + sample + '.bed; \
                      /mnt/lustre/home/cusanovich/Programs/BEDTools/bin/coverageBed -a ' + sample + '.bed -b /mnt/lustre/data/users/cusanovich/References/hg19ProteinCodingEnsemblExonsMergedNonoverlapping.bed > ' + sample + '.exoncounts.txt; \
                      python /mnt/lustre/data/users/cusanovich/RNAseq_scripts/exoncombiner.py ' + sample + '.exoncounts.txt ' + sample + '.genecounts.txt; \
                      rm ' + sample + '.quality.merged.bam; \
                      rm ' + sample + '.junction.quality.merged.bam; \
                      rm ' + sample + '.junction.quality.merged.sam; \
                      rm ' + sample + '.junction.bed" | qsub -l h_vmem=8g -o ' + outdir + fcid + ' -e ' + outdir + fcid + ' -N "'+ sample.split('/')[-1] + '.mergers"'
#            print merger
            merging = subprocess.Popen(merger, shell=True)
            merging.communicate()
            x += 1


