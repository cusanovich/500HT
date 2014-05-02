import os
import genecounter
import pysam
import numpy
import sys
file_extension = 'quality.sort.bam'
job_id = 'genotyper'

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
#fcno = '3'
#fcid = 'FlowCell' + fcno + '/'
sample = sys.argv[1]
fcid = sample.split('/')[-3]
currsample = sample.split('/')[-2]
currlane = 'fc' + fcno + '.' + currsample + '.' + sample.split('/')[-1]
currright = 0
currwrong = 0
currtotal = 0
currsnps = 0
currind = findivs.index(currsample)
print currind
reversed = dict({'A':'T', 'T':'A', 'G':'C', 'C':'G'})
genos = open("/mnt/lustre/home/cusanovich/oldhome_nb/Hutterite_Heritability/gemma/Genotypes/hutt.3chip.500HT.tped","r")
homogeno = {}
for line in genos:
    if line == '':
        continue
    line = line.strip('\n').split(' ')
    Aind = currind*2 + 4
    Bind = currind*2 + 5
#    Aind = currind*2 + 6  #This is to pull the wrong alleles, just to check my null
#    Bind = currind*2 + 7  #This is to pull the wrong alleles, just to check my null
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

print len(homogeno.keys())
skipper = 0
inder = sample.split('/')[-1] + '.quality.sort.bam.bai'
filer = os.listdir('/mnt/lustre/data/users/cusanovich/500HTRNAseq/' + fcid + str(findivs[currind]) + '/')
for file in filer:
    if file == inder:
        skipper = 1

if skipper != 1:
    print 'Indexing...'
    pysam.index(sample + '.quality.sort.bam')

samfile = pysam.Samfile(sample + '.quality.sort.bam','rb')
x = 0
for geno in homogeno.keys():
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
#            print geno + ' has ' + str(ner) + ' reads mapped.'
            currtotal = currtotal + ner
            currsnps += 1
            for pileupread in pileupcolumn.pileups:
                seqer.append(pileupread.alignment.seq[pileupread.qpos])
    for read in seqer:
        if read == homogeno[geno][0] or read == reversed[homogeno[geno][0]]:
            currright += 1
        if read == homogeno[geno][1] or read == reversed[homogeno[geno][1]]:
            currwrong += 1
    if x%10000 == 0:
        print str(x) + ' SNPs analyzed...'
    x += 1

samfile.close()
print currlane + '\t' + str(currright) + '\t' + str(currwrong) + '\t' + str(currtotal) + '\t' + str(currsnps) + '\t' + str(float(currright)/float(currtotal))
oldind = currind

outfile = open(sample + '.genotype.report.txt','w')
print >> outfile, str(currright) + '\n' + str(currwrong) + '\n' + str(currtotal) + '\n' + str(currsnps) + '\n' + str(float(currright)/float(currtotal))
outfile.close()