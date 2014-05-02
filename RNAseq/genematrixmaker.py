#import sys
#import os
import subprocess
import genecounter
import numpy

if(__name__=='__main__'):
    x = 0
    y = 0
    outdir = '/mnt/lustre/data/users/cusanovich/RNAseq_scripts/bedmakerOE/'
    header = ['Chr', 'Start', 'End', 'Gene', 'ExonicLength']
    for j in range(1,12):
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
            sampler = open(sample + '.gccor.genecounts.txt','r')
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
                    length.append(line[5])
                counts.append(line[4])
            chrm = numpy.array(chrm)
            #print counts[0:10]
            #print gene[0:10]
            if x == 0:
                genecounting = numpy.column_stack((numpy.array(chrm)[:,numpy.newaxis],numpy.array(start)[:,numpy.newaxis],
                                                   numpy.array(end)[:,numpy.newaxis],numpy.array(gene)[:,numpy.newaxis],
                                                   numpy.array(length)[:,numpy.newaxis],numpy.array(counts)[:,numpy.newaxis],))
            else:
                genecounting = numpy.column_stack((genecounting,numpy.array(counts)[:,numpy.newaxis]))
            x = 1
    genecounting = numpy.row_stack((header,genecounting))
    numpy.savetxt('gccor.genecount.matrix', genecounting, fmt = '%s', delimiter = '\t', newline = '\n')
    #numpy.savetxt('fc' + fcno + 'genecount.matrix', genecounting, fmt = '%s', delimiter = '\t', newline = '\n')
