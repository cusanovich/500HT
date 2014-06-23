#!/usr/bin/env python
import subprocess

print 'Converting GTF to BED...'
cager = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/wgEncodeRikenCageGm12878CellPapTssGencV7.gtf','r')
outter = open('/mnt/lustre/home/cusanovich/500HT/ensemblCAGETSS.bed','w')
for line in cager:
    liner = line.strip().split()
    chrom = liner[0]
    start = str(int(liner[3]) - 1)
    end = liner[4]
    score = liner[5]
    strand = liner[6]
    gene = liner[9].strip('"').split('.')[0]
    biotype = liner[13].strip('"').split(',')[0]
    print >> outter, chrom + '\t' + start + '\t' + end + '\t' + gene + '\t' + score + '\t' + strand + '\t' + biotype

cager.close()
outter.close()

print 'Collecting expressed genes...'
exprgenes = open('/mnt/lustre/home/cusanovich/500HT/genenames.500ht.txt','r')
exprs = exprgenes.readlines()
exprs = [x.strip() for x in exprs]

print 'Consolidating genes...'
bed = open('/mnt/lustre/home/cusanovich/500HT/ensemblCAGETSS.bed','r')
newbedfile = '/mnt/lustre/home/cusanovich/500HT/ensemblCAGETSS_RNAseq.bed'
newbed = open(newbedfile,'w')
sortbed = '/mnt/lustre/home/cusanovich/500HT/ensemblCAGETSS_RNAseq_sorted.bed'
print 'Building gene dictionary...'
chrs = {}
starts = {}
ends = {}
scores = {}
strands = {}
biotypes = {}

for line in bed:
    if line == '':
        continue
    liner = line.strip().split()
    try:
        chrs[liner[3]].append(liner[0])
    except KeyError:
        chrs[liner[3]] = [liner[0]]
    try:
        starts[liner[3]].append(int(liner[1]))
    except KeyError:
        starts[liner[3]] = [int(liner[1])]
    try:
        ends[liner[3]].append(int(liner[2]))
    except KeyError:
        ends[liner[3]] = [int(liner[2])]
    try:
        scores[liner[3]].append(float(liner[4]))
    except KeyError:
        scores[liner[3]] = [float(liner[4])]
    try:
        strands[liner[3]].append(liner[5])
    except KeyError:
        strands[liner[3]] = [liner[5]]
    try:
        biotypes[liner[3]].append(liner[6])
    except KeyError:
        biotypes[liner[3]] = [liner[6]]

bed.close()

print 'Dictionary contains ' + str(len(chrs.keys())) + ' genes...'
print 'Collapsing records...'
j=0
k=0
m=0
for bedrecord in chrs.keys():
    if bedrecord not in exprs:
        continue
    if len(chrs[bedrecord]) == 1 and biotypes[bedrecord] == 'protein_coding':
        print >> newbed, chrs[bedrecord][0] + '\t' + str(starts[bedrecord][0]) + '\t' + str(ends[bedrecord][0]) + '\t' + bedrecord + '\t' + str(scores[bedrecord][0]) + '\t' + strands[bedrecord][0]
        j += 1
        continue
    chra = chrs[bedrecord][0]
    stranda = strands[bedrecord][0]
    if biotypes[bedrecord].count('protein_coding') == 0:
        winning = max(scores[bedrecord])
        winners = []
        count = -1
        for score in scores[bedrecord]:
            count += 1
            if score == winning:
                winners.append(count)
        newstarts = []
        wins = [str(winning)]
        for g in range(len(winners)):
            newstarts.append(starts[bedrecord][winners[g]])
            if len(winners) > 1:
                wins.append(str(starts[bedrecord][winners[g]]))
        newend = int(min(newstarts))
        if stranda == "-":
            newend = int(max(newstarts))
        newstart = newend - 1
        newscore = "_".join(wins)
        #print >> newbed, chra + '\t' + str(newstart) + '\t' + str(newend) + '\t' + bedrecord + '\t' + newscore + '\t' + stranda
        k += 1
        continue
    coding = []
    count = -1
    for bios in biotypes[bedrecord]:
        count += 1
        if bios == 'protein_coding':
            coding.append(count)
    newstarts = []
    newscores = []
    for g in range(len(coding)):
    	newstarts.append(starts[bedrecord][coding[g]])
    	newscores.append(scores[bedrecord][coding[g]])
    winning = max(newscores)
    winners = []
    count = -1
    for score in newscores:
        count += 1
        if score == winning:
            winners.append(count)
    newerstarts = []
    wins = [str(winning)]
    for g in range(len(winners)):
        newerstarts.append(newstarts[winners[g]])
        if len(winners) > 1:
            wins.append(str(newstarts[winners[g]]))
    newend = int(min(newstarts))
    if stranda == "-":
        newend = int(max(newstarts))
    newstart = newend - 1
    newscore = "_".join(wins)
    print >> newbed, chra + '\t' + str(newstart) + '\t' + str(newend) + '\t' + bedrecord + '\t' + newscore + '\t' + stranda
    m += 1


newbed.close()
print 'There were ' + str(j) + ' genes with 1 TSS.'
print 'There were ' + str(k) + ' genes with no protein coding TSS.'
print 'There were ' + str(m) + ' genes with at least 1 protein coding TSS.'

print 'Sorting records...'
sorter = 'sort -k1,1 -k2,2n -k3,3n -k4,4 ' + newbedfile + ' > ' + sortbed
sorting = subprocess.Popen(sorter,shell=True)
sorting.wait()
