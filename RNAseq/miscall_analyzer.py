#!/usr/bin/env python
import os
import pysam
import numpy
import sys
file_extension = 'bed'
job_id = 'miscalls'
gb_request = '16g'

findivs = open("/mnt/lustre/home/cusanovich/oldhome_nb/Hutterite_Heritability/gemma/Genotypes/500HT.geno.findivs.txt","r").read().splitlines()
alleles = {}
allelic = open("/mnt/lustre/home/cusanovich/oldhome_nb/Hutterite_Heritability/gemma/Genotypes/hutt.3chip.500HT.bim","r")
for line in allelic:
    if line == '':
        continue
    line = line.strip('\n').split('\t')
    alleles[line[1]] = [line[0],line[3],line[4],line[5]]

allelic.close()
oldind = 'none'
sample = sys.argv[1]
currsample = sample.split('/')[-2]
try:
    currind = findivs.index(currsample)
except ValueError:
    sys.stderr.write(currsample + ' is missing!')
    sys.exit(1)
concordance = open(sample + ".genotype.report.txt","r").read().splitlines()[4]
#if float(concordance) < 0.9 or float(concordance) > 0.97:
#    sys.exit(1)
fcid = sample.split('/')[-3]
currind = findivs.index(currsample)
reversed = dict({'A':'T', 'T':'A', 'G':'C', 'C':'G'})
genos = open("/mnt/lustre/home/cusanovich/oldhome_nb/Hutterite_Heritability/gemma/Genotypes/hutt.3chip.500HT.tped","r")
homogeno = {}
for line in genos:
    if line == '':
        continue
    line = line.strip('\n').split(' ')
    Aind = currind*2 + 4
    Bind = currind*2 + 5
    alleleA = line[Aind]
    alleleB = line[Bind]
    if alleleA != alleleB:
        continue
    if alleleA == '0':
        continue
    if alleleA == alleles[line[1]][2]:
        wrongallele = alleles[line[1]][3]
    else:
        wrongallele = alleles[line[1]][2]
    if alleleA == reversed[wrongallele]:
        continue
    homogeno[line[1]] = [alleleA,wrongallele]

skipper = 0
inder = sample.split('/')[-1] + '.quality.sort.bam.bai'
filer = os.listdir('/mnt/lustre/data/users/cusanovich/500HTRNAseq/' + fcid + '/' + str(findivs[currind]) + '/')
for file in filer:
    if file == inder:
        skipper = 1

if skipper != 1:
    sys.stderr.write('How come this .bam is not indexed?')
    sys.exit(1)
samfile = pysam.Samfile(sample + '.quality.sort.bam','rb')
outfile = open(sample + '.genotype.bysnp.report.txt','w')
print >> outfile, 'rs\tNoRight\tNoWrong\tTotal'
x = 0
for geno in homogeno.keys():
    currright = 0
    currwrong = 0
    seqer = []
    ner = 0
    currchr = 'chr' + alleles[geno][0]
    if currchr == 'chr23':
        currchr = 'chrX'
    currstart = int(alleles[geno][1])-1
    currend = int(alleles[geno][1])
    for pileupcolumn in samfile.pileup(currchr, currstart, currend):
        if pileupcolumn.pos == currstart:
            ner = pileupcolumn.n
            for pileupread in pileupcolumn.pileups:
                seqer.append(pileupread.alignment.seq[pileupread.qpos])
    for read in seqer:
        if read == homogeno[geno][0] or read == reversed[homogeno[geno][0]]:
            currright += 1
        if read == homogeno[geno][1] or read == reversed[homogeno[geno][1]]:
            currwrong += 1
    if ner > 0:
        print >> outfile, geno + '\t' + str(currright) + '\t' + str(currwrong) + '\t' + str(ner)
    if x%10000 == 0:
        sys.stdout.write(str(x) + ' SNPs analyzed...\n')
    if x%1000 == 0:
        sys.stdout.write('.')
    x += 1

samfile.close()
outfile.close()
