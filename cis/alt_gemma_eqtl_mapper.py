#!/usr/bin/env python
####Import needed modules
import os
import sys
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
import subprocess
import numpy
import random
import gzip
import time

####Check that proper arguments are supplied
#if len(sys.argv) != 3:
#	print "Usage = python gemma_eqtl_mapper.py chr pc"
#	sys.exit()

chrm = sys.argv[1]
pcs = sys.argv[2]
#chrm = 'chr22'
#pcs = 0
genodir = '/mnt/lustre/home/cusanovich/500HT/3chip/'
hmdir = '/mnt/lustre/home/cusanovich/'
currfiles = genodir + "curr_" + chrm + "_pc" + str(pcs)

####Define handy functions: 'ifier' helps to communicate with the shell from python
#'matrix_reader' reads in large datasets more quickly and efficiently
def ifier(commander):
	ify = subprocess.Popen(commander,shell=True)
	ify.wait()

def matrix_reader(matrix_file,sep="\t",dtype='|S20'):
	linecounter = subprocess.Popen('wc -l ' + matrix_file, shell=True, stdout=subprocess.PIPE)
	linecount = int(linecounter.communicate()[0].strip().split()[0])
	columncounter = subprocess.Popen('awk -F"' + sep + '" \'{print NF;exit}\' ' + matrix_file, shell=True, stdout=subprocess.PIPE)
	columncount = int(columncounter.communicate()[0].strip().split()[0])
	raws = numpy.zeros((linecount,columncount),dtype=dtype)
	rawin = open(matrix_file,'r')
	for i,line in enumerate(rawin):
		raws[i,:] = line.strip().split()
	rawin.close()
	return raws

####Added this because Plink seems to crash if multiple jobs are trying to access the same plink file at the same time
time.sleep(random.randint(1,150))
####Read in map of which columns of expression matrix have which SNPs in cis
mastercols = matrix_reader('/mnt/lustre/home/cusanovich/500HT/hutt.3chip.mastercols.txt',dtype='|S15')
####Build dictionaries (like a hash in perl) to (a) keep track of which SNPs go with which genes,
####(b) keep track of which column of the expression matrix belongs to each gene, and (c) which chromosome each gene is on
masterdic = {}
exprcoldic = {}
chrmdic = {}
for i in range(mastercols.shape[0]):
	try:
		masterdic[mastercols[i,0]].append(mastercols[i,1])
	except KeyError:
		masterdic[mastercols[i,0]] = [mastercols[i,1]]
		exprcoldic[mastercols[i,0]] = mastercols[i,2]
		chrmdic[mastercols[i,0]] = mastercols[i,4]

pvals = []
genes = []
snps = []
winnerdic = {}
pcdic = {}
#completedgenes = 0
chr22ers = []
for gene in masterdic.keys():
	if chrmdic[gene] == chrm:
		chr22ers.append(gene)

#t0 = time.time()

####Loop through each gene to calculate eQTL p-values
for gene in masterdic.keys():
#for gene in chr22ers[0:1]:
#for gene in masterdic.keys()[1:500]:
	if chrmdic[gene] != chrm:
		continue
#	if completedgenes == 20:
#		continue
	print gene
	#print "Extracting SNPs..."
	####Pull genotypes for all cis SNPs for the current gene
	extracter = "plink --noweb --bfile " + genodir + "hutt.3chip --from " + masterdic[gene][0] + " --to " + masterdic[gene][-1] + " --make-bed --out " + currfiles
	ifier(extracter)
	#if findivs[0] not in open(currfiles + '.fam','r').readline():
	#	print "Warning! Genotype findivs do not match phenotype findivs!"
	#	sys.exit()
	#print "Prepping GEMMA inputs..."
	####Creating files to feed into GEMMA
	famer = ('cut -f1-5 -d" " ' + currfiles + '.fam > ' + currfiles + '.header.fam; cut -f' + str(int(exprcoldic[gene]) + 6) + ' -d" " ' + genodir +
		'hutt.3chip.fam.update > ' + currfiles + '.tailer.fam; paste -d" " ' + currfiles + '.header.fam ' + currfiles + '.tailer.fam > ' +
		currfiles + '.fam; rm ' + currfiles + '.header.fam; rm ' + currfiles + '.tailer.fam')
	ifier(famer)
	#print "Running GEMMA..."
	####Calling GEMMA with the current cis SNPs and the current gene
	gemmer = ('cd ' + genodir + '; ' + hmdir + 'Programs/gemma0.94 -bfile curr_' + chrm + '_pc' + str(pcs) + ' -km 2 -k ' + hmdir + '500HT/addSNP.500ht.txt -c ' + hmdir + '500HT/Exprs/qqnorm.500ht.3chip_order.pc' + str(pcs) + ' -lmm 2 -o curr_' + chrm + '_pc' + str(pcs))
	ifier(gemmer)
	currresults = open(genodir + '/output/curr_' + chrm + '_pc' + str(pcs) + '.assoc.txt','r')
	####Collecting best cis SNP from GEMMA results
	pmin = 1.1
	snpmin = ''
	for line in currresults:
		if 'p_lrt' in line:
			continue
		liner = line.strip().split()
		if liner[5] == 'nan':
			continue
		genes.append(gene)
		snps.append(liner[1])
		pvals.append(float(liner[5]))
		if float(liner[5]) < pmin:
			pmin = float(liner[5])
			snpmin = liner[1]
	#print "Starting up permutations..."
	winnerperms = 0
	broken = 0
	findic = {}
	####Setting up for permutations
	for perm in range(10000):
		if broken == 1:
			continue
		####Building dictionaries to keep track of which expression and pcs belong to which individual
		if perm == 0:
			try:
				test = findic.keys()[0]
			except IndexError:
				finder = open(currfiles + '.fam','r')
				findivs = []
				randfindivs = []
				for line in finder:
					liner = line.strip().split()
					findic[liner[1]] = liner
					findivs.append(liner[1])
					randfindivs.append(liner[1])
				finder.close()
			try:
				test = pcdic.keys()[0]
			except IndexError:
				pcfile = open(hmdir + '500HT/Exprs/qqnorm.500ht.3chip_order.pc' + str(pcs),'r')
				for z,line in enumerate(pcfile):
					pcdic[findivs[z]] = line.strip().split()
				pcfile.close()
		####Randomize findivs and ask plink to rearrange the genotype files
		random.shuffle(randfindivs)
		rearrange = open(currfiles + '.findivs.rearrange.txt','w')
		for z in range(len(findivs)):
			print >> rearrange, 'HUTTERITES\t' + findivs[z] + '\tHUTTERITES\t' + randfindivs[z]
		rearrange.close()
		plinker = 'plink --noweb --bfile ' + currfiles + ' --update-ids ' + currfiles + '.findivs.rearrange.txt --make-bed --out ' + genodir + 'perm_curr_' + chrm + '_pc' + str(pcs)
		ifier(plinker)
		####Put the expression and PC data back in order so that the link between genotypes and phenotypes is broken
		permfam = open(genodir + 'perm_curr_' + chrm + '_pc' + str(pcs) + '.fam','w')
		permpc = open(genodir + 'perm_curr_' + chrm + '_pc' + str(pcs) + '.pcs','w')
		for z in range(len(findivs)):
			print >> permfam, " ".join(findic[randfindivs[z]])
			print >> permpc, " ".join(pcdic[randfindivs[z]])
		permfam.close()
		permpc.close()
		print str(winnerperms)
		####Run GEMMA on permuted data
		gemmer = ('cd ' + genodir + '; ' + hmdir + 'Programs/gemma0.94 -bfile perm_curr_' + chrm + '_pc' + str(pcs) + ' -km 2 -k ' + hmdir + '500HT/addSNP.500ht.txt -c ' + genodir + 'perm_curr_' + chrm + '_pc' + str(pcs) + '.pcs' + ' -lmm 2 -o perm_curr_' + chrm + '_pc' + str(pcs))
		ifier(gemmer)
		####Count if any cis permuted p-value is better than observed cis p-value
		perms = []
		permer = open(genodir + 'output/perm_curr_' + chrm + '_pc' + str(pcs) + '.assoc.txt','r')
		permlow = 1
		for line in permer:
			liner = line.strip().split()
			if liner[5] == 'nan' or 'p_lrt' in liner[5]:
				continue
			if float(liner[5]) < permlow:
				permlow = float(liner[5])
		permer.close()
		if permlow <= pmin:
			winnerperms += 1
		if winnerperms == 10:
			genep = 11/random.uniform(perm+2,perm+3)
			broken = 1
	####Calculate permuted p-value
	if broken == 0:
		genep = (winnerperms + 1)/float(10001)
	####Record permuted p-value
	winnerdic[gene] = [snpmin, pmin, genep]
	cleanup = 'rm ' + hmdir + '500HT/3chip/*curr_' + chrm + '_pc' + str(pcs) + '*'
	ifier(cleanup)
#	completedgenes += 1

#t1 = time.time()
#print t1-t0

print "Writing results..."
####Write out all GEMMA p-values
aller = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.3chip.gemma.eqtls.txt','w')
for i in range(len(genes)):
	print >> aller, '{0}\t{1}\t{2:.4g}'.format(genes[i],snps[i],pvals[i])

aller.close()

####Write out best SNP and permuted p-value for each gene
winners = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.3chip.gemma.chosen.txt','w')
for gene in sorted(winnerdic.keys()):
	print >> winners, '{0}\t{1[0]}\t{1[1]:.4g}\t{1[2]:.4g}'.format(gene,winnerdic[gene])

winners.close()

doners = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.done','w')
doners.close()
