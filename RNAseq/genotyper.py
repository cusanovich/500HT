import os
import genecounter
import pysam
import numpy

if(__name__=='__main__'):
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
    reversing = dict({'A':'T', 'T':'A', 'G':'C', 'C':'G'})
#    right = []
#    wrong = []
#    total = []
#    results = {}
#    reruns = [7,9]
#    for j in reruns:
    for j in range(11,12):
        fcno = str(j)
        fcid = 'FlowCell' + fcno + '/'
        print fcid + '...\n...\n...\n...\n...\n...\n...\n'
        finlist = genecounter.findiv_grabber(fcid)
#        print finlist[0]

        samples = []
        for id in finlist:
            countslist = genecounter.lane_grabber(fcid,id)
            samples.extend(countslist)
#        print samples[0:11]

        for sample in samples:
#            print sample
            if os.path.isfile(sample + '.genotype.report.txt'):
                continue
            currsample = sample.split('/')[-2]
            currlane = 'fc' + fcno + '.' + currsample + '.' + sample.split('/')[-1]
            print currlane
#            print currsample
            currright = 0
            currwrong = 0
            currtotal = 0
            currsnps = 0
            try:
                currind = findivs.index(currsample)
                print currind
            except ValueError:
                print currsample + ' is missing!'
                continue
            if oldind != currind:
                print 'Collecting SNPs...'
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
                    if alleleA == reversing[alleleB]:
                        continue
                    if alleleA == alleles[line[1]][2]:
                        wrongallele = alleles[line[1]][3]
                    else:
                        wrongallele = alleles[line[1]][2]
                    if alleleA == reversing[wrongallele]:
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
            print 'Opening SAM file...'
            samfile = pysam.Samfile(sample + '.quality.sort.bam','rb')
            print 'Collecting genotype calls...'
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
                        #print geno + ' has ' + str(ner) + ' reads.'
                        currtotal = currtotal + ner
                        currsnps += 1
                        for pileupread in pileupcolumn.pileups:
                            seqer.append(pileupread.alignment.seq[pileupread.qpos])
                            #currtotal +=1
                for read in seqer:
                    if read == homogeno[geno][0] or read == reversing[homogeno[geno][0]]:
                        currright += 1
                    if read == homogeno[geno][1] or read == reversing[homogeno[geno][1]]:
                        currwrong += 1
            samfile.close()
            print currlane + '\t' + str(currright) + '\t' + str(currwrong) + '\t' + str(currtotal) + '\t' + str(currsnps) + '\t' + str(float(currright)/float(currtotal))
#            results[currlane] = [currright,currwrong,currtotal,currsnps,float(currright)/float(currtotal)]
            oldind = currind
            outfile = open(sample + '.genotype.report.txt','w')
            print >> outfile, str(currright) + '\n' + str(currwrong) + '\n' + str(currtotal) + '\n' + str(currsnps) + '\n' + str(float(currright)/float(currtotal))
#            print >> outfile, "Lane\tMatching\tWrong\tTotal\tNoSNPs\tPercent"
#            for key in sorted(results.iterkeys()):
#                print >> outfile, key + '\t' + str(results[key][0]) + '\t' + str(results[key][1]) + '\t' + str(results[key][2]) + '\t' + str(results[key][3]) + '\t' + str(results[key][4])
            outfile.close()
            print 'Done with sample.'