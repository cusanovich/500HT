import subprocess
import genecounter
import numpy

# 1: get all count files
# 2: read in records for each individual
# 3: make a matrix with a row for each type of mapping and a column for each individual

def digitizer(number):
    number = str(number)
    if len(number) == 1:
        number = '0' + number
    return number

if(__name__=='__main__'):
    rownames = ['ReadSource', 'Total', 'Demultied','Saved', 'Mapped', 'JuncMapped', 'SaveMapped', 'JuncSaved']
    reads = {}
    badlist = []
    qcfile = open('/mnt/lustre/home/cusanovich/data/RNAseq_scripts/genotypes.matrix.update')
    for line in qcfile:
        if 'SampleID' in line:
            continue
        if float(line.strip().split()[4]) < 0.975:
            badlist.append(line.strip().split()[0])
    print 'There are ' + str(len(badlist)) + ' bad samples (by genotype).'
    for j in range(1,12):
        fcno = str(j)
        fcid = 'FlowCell' + fcno + '/'
        print fcid
        finlist = genecounter.findiv_grabber(fcid)
        print finlist[0]

        samples = []
        for id in finlist:#['163802','168851']:#
            countslist = genecounter.lane_grabber(fcid,id)
            for countfile in countslist:
                countfields = countfile.split('/')
                sampler = countfields[8]
                fc = digitizer(countfields[7].split('FlowCell')[1])
                laner = digitizer(countfields[9].split('.')[0].split('_')[1])
                inder = digitizer(countfields[9].split('.')[1].split('_')[1])
                sampleid = sampler + '.' + fc + laner + inder
                if sampler == '160462':
                    print sampleid
                if sampleid in badlist:
                    continue
                samples.append(countfile)
#        print samples[0:11]
        print 'There are ' + str(len(samples)) + ' good samples.'
        for sample in samples:
            samp = sample.split('/')[-2]
            print samp
            sampler = open(sample + '.counts.txt','r')
            counts = sampler.readlines()
            tot = int(counts[0].strip().split('\t')[1]) + int(counts[5].strip().split('\t')[1])
            totd = int(counts[0].strip().split('\t')[1])
            tots = int(counts[5].strip().split('\t')[1])
            mapped = int(counts[1].strip().split('\t')[1])
            juncmapped = int(counts[3].strip().split('\t')[1])
            savemapped = int(counts[6].strip().split('\t')[1])
            juncsaved = int(counts[7].strip().split('\t')[1])
            if samp not in reads.keys():
                reads[samp] = [tot,totd,tots,mapped,juncmapped,savemapped,juncsaved]
            else:
                reads[samp][0] += tot
                reads[samp][1] += totd
                reads[samp][2] += tots
                reads[samp][3] += mapped
                reads[samp][4] += juncmapped
                reads[samp][5] += savemapped
                reads[samp][6] += juncsaved
#            print reads

    print str(len(reads.keys())) + ' samples have good reads.'
    readings = numpy.array(rownames)[:,numpy.newaxis]
#    print readings
    for key in sorted(reads.iterkeys()):
        currlist = [str(key), str(reads[key][0]), str(reads[key][1]), str(reads[key][2]), str(reads[key][3]), str(reads[key][4]), str(reads[key][5]), str(reads[key][6])]
#        print currlist
#        print readings.shape
        readings = numpy.column_stack((readings,currlist))
        numpy.savetxt('reads.matrix', readings, fmt = '%s', delimiter = '\t', newline = '\n')


# Example File Format
# Total:	15372199
# Mapped:	10872866
# Unmapped:	2650833
# JunctionMapped:	1559469
# JunctionUnmapped:	1016078
# Saved:	494853
# SavedMapped:	337104
# SavedUnmapped:	100246
# SavedJunctionMapped:	49213
# SavedJunctionUnmapped:	48491