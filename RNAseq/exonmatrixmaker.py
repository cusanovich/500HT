#import sys
#import os
import sys
import subprocess
import genecounter
import numpy
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
from DarrenTools import ifier

if(__name__=='__main__'):
    x = 0
    y = 0
    z = 1
    fclist = []
    outdir = '/mnt/lustre/data/users/cusanovich/RNAseq_scripts/bedmakerOE/'
    header = ['Chr', 'Start', 'End', 'Gene', 'ExonicLength']
    for j in range(1,12):
        if x == 1:
            header = []
            z = 0
        fcno = str(j)
        fcid = 'FlowCell' + fcno + '/'
        print fcid
        finlist = genecounter.findiv_grabber(fcid)
        print finlist[0]

        samples = []
        for ids in finlist:#['163802','168851']:#
            countslist = genecounter.lane_grabber(fcid,ids)
            samples.extend(countslist)
#        print samples[0:11]

        for sample in samples:
            #print sample
            fcnum = '.fc0' + str(fcno)
            if int(fcno) > 9:
                fcnum = '.fc' + str(fcno)
            samp = sample.split('/')[-2] + '.' + sample.split('/')[-1].split('.')[0] + fcnum
            header.append(samp)
            print samp
            if y == 0:
                chrm = []
                start = []
                end = []
                gene = []
                length = []
                y = 1
            counts = []
            try:
                sampler = open(sample + '.exoncounts.txt','r')
#            sampler = open(sample + '.genecounts.txt','r')
                for line in sampler:
                    if line == '':
                        continue
                    line = line.strip('\n').split('\t')
                    if line[0] == 'Chr':
                        continue
                    if x == 0:
                        chrm.append(line[0])
                        start.append(line[1])
                        end.append(line[2])
                        gene.append(line[3])
                        length.append(line[6])
                    counts.append(line[4])
                chrm = numpy.array(chrm)
                if x == 0:
                    genecounting = numpy.column_stack((numpy.array(chrm)[:,numpy.newaxis],numpy.array(start)[:,numpy.newaxis],
                                                       numpy.array(end)[:,numpy.newaxis],numpy.array(gene)[:,numpy.newaxis],
                                                       numpy.array(length)[:,numpy.newaxis],numpy.array(counts)[:,numpy.newaxis],))
                else:
                    if z == 0:
                        genecounting = numpy.array(counts)[:,numpy.newaxis]
                        z = 1
                    else:
                        genecounting = numpy.column_stack((genecounting,numpy.array(counts)[:,numpy.newaxis]))
                x = 1
            #Code to remove corrupted exoncounts files and regenerate...
            except IOError:
                print 'Fixing ' + samp + '...'
                ifier('rm ' + sample + '.exoncounts.txt')
                cleanup = '/mnt/lustre/home/cusanovich/Programs/samtools/samtools merge ' + sample + '.quality.merged.bam ' + sample + '.quality.sort.bam ' + sample + '.saved.quality.sort.bam; \
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
                           rm ' + sample + '.junction.bed; \
                           python 500HT/Scripts/RNAseq/exonmatrixmaker.py;'
                ifier(cleanup)
                print samp + ' fixed.'
                for line in sampler:
                    if line == '':
                        continue
                    line = line.strip('\n').split('\t')
                    if line[0] == 'Chr':
                        continue
                    if x == 0:
                        chrm.append(line[0])
                        start.append(line[1])
                        end.append(line[2])
                        gene.append(line[3])
                        length.append(line[6])
                    counts.append(line[4])
                chrm = numpy.array(chrm)
            #print counts[0:10]
            #print gene[0:10]
                if x == 0:
                    genecounting = numpy.column_stack((numpy.array(chrm)[:,numpy.newaxis],numpy.array(start)[:,numpy.newaxis],
                                                       numpy.array(end)[:,numpy.newaxis],numpy.array(gene)[:,numpy.newaxis],
                                                       numpy.array(length)[:,numpy.newaxis],numpy.array(counts)[:,numpy.newaxis],))
                else:
                    if z == 0:
                        genecounting = numpy.array(counts)[:,numpy.newaxis]
                        z = 1
                    else:
                        genecounting = numpy.column_stack((genecounting,numpy.array(counts)[:,numpy.newaxis]))
                x = 1
        print len(header)
        print genecounting.shape
        print header[0:10]
        genecounting = numpy.row_stack((header,genecounting))
        numpy.savetxt('fc' + fcno + '.exoncount.matrix', genecounting, fmt = '%s', delimiter = '\t', newline = '\n')
        fclist.append('fc' + fcno + '.exoncount.matrix')
    print "Joining flow cells..."
    combo = 'paste ' + " ".join(fclist) + ' > master.exoncount.matrix'
    ifier(combo)
    #genecounting = numpy.row_stack((header,genecounting))
    #numpy.savetxt('exoncount.matrix', genecounting, fmt = '%s', delimiter = '\t', newline = '\n')
    #numpy.savetxt('fc' + fcno + 'genecount.matrix', genecounting, fmt = '%s', delimiter = '\t', newline = '\n')
