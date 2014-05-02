import subprocess
import genecounter
import numpy
import os.path

# 1: get all genotype count files
# 2: read in records for each sample
# 3: make a matrix with a row for each type of mapping and a column for each sample
# This will name samples with the following naming convention:
#     findiv.(two digit fcno)(two digit laneID)(two digit IndexID)
# For example, 30581.070612 would be sample 30581 run on flow cell 7, lane 6, tagged with index 12.

if(__name__=='__main__'):
    rownames = ['SampleID', 'RightReads', 'WrongReads','TotalReads','PercentRight','No.SNPs']
    reads = {}
    for j in range(1,12):
        fcno = str(j)
        fcid = 'FlowCell' + fcno + '/'
        print fcid
        finlist = genecounter.findiv_grabber(fcid)
        print finlist[0]

        samples = []
        for id in finlist:#['163802','168851']:#
            countslist = genecounter.lane_grabber(fcid,id)
            samples.extend(countslist)
#        print samples[0:11]

        for sample in samples:
            samp = str(sample.split('/')[-2])
            laner = str(sample.split('/')[-1].split('.')[0].split('_')[1])
            indexer = str(sample.split('/')[-1].split('.')[1].split('_')[1])
#            print samp
            if len(fcno) < 2:
                fcsamp = '0' + fcno
            else:
                fcsamp = fcno
            if len(laner) < 2:
                laner = '0' + laner
            if len(indexer) < 2:
                indexer = '0' + indexer
            fullsamp = samp + '.' + fcsamp + laner + indexer
#            print fullsamp
            try:
                sampler = open(sample + '.genotype.report.txt','r')
                counts = sampler.readlines()
                right = str(counts[0].strip())
                wrong = str(counts[1].strip())
                tot = str(counts[2].strip())
                snps = str(counts[3].strip())
                perc = str(counts[4].strip())
                reads[fullsamp] = [right,wrong,tot,perc,snps]
            except IOError:
                print sample + ' is not genotyped!'

#    readings = numpy.array(rownames)[:,numpy.newaxis]
    readings = numpy.array(rownames)
    for key in sorted(reads.iterkeys()):
        currlist = [str(key), str(reads[key][0]), str(reads[key][1]), str(reads[key][2]), str(reads[key][3]), str(reads[key][4])]
#        print currlist
#        print readings.shape
        readings = numpy.row_stack((readings,currlist))
        numpy.savetxt('genotypes.matrix', readings, fmt = '%s', delimiter = '\t', newline = '\n')


# Example File Format
# 20603		#Total No. reads matching the correct allele or its reverse complement
# 296			#Total No. reads matching the wrong allele or its reverse complement
# 20899		#Total No. reads overlapping a testable SNP
# 3999			#Total No. SNPs covered and testable
# 0.985836642902	#Total No. correct reads / total No. reads (i.e. line 1/line 3)