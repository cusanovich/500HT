#!/usr/bin/env python
import os
import sys
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
import subprocess
import numpy
import random
import gzip
import time

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
mapper = '3chip'
distance = '150kb'

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

print "Loading expression..."
mastercols = matrix_reader('/mnt/lustre/home/cusanovich/500HT/hutt.' + mapper + '.' + distance + '.mastercols.txt',dtype='|S15')
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

####Build a dictionary to reference the genomic coordinates of each SNP
print "Loading SNP annotations..."
snpdic = {}
snpbed = open('/mnt/lustre/home/cusanovich/500HT/hutt.3chip.hg19.bed','r')
for line in snpbed:
	liner = line.strip().split()
	snpdic[liner[3]] = liner[0:3]

cover = ('cp ' + hmdir + '500HT/addSNP.500ht.' + mapper + '_order.square.txt ' + currfiles + '.square.txt; cp ' + hmdir + '500HT/Exprs/qqnorm.500ht.' + mapper + '_order.pc' + str(pcs)) + ' ' + currfiles + '.pcs.txt'
ifier(cover)


#completedgenes = 0
#chr22ers = []
#for gene in masterdic.keys():
#	if chrmdic[gene] == chrm:
#		chr22ers.append(gene)


pvals = []
genes = []
snps = []
winnerdic = {}
#t0 = time.time()

for gene in masterdic.keys():
	if chrmdic[gene] != chrm:
		continue
#	if completedgenes == 20:
#		continue
	print gene
	currbimbam = open(currfiles + '.bimbam','r')
	####Pull genotypes for the SNPs in cis, if genotypes not already in dictionary: go to geno file and pull in appropriate data
	for i, snp in enumerate(masterdic[gene]):
		try:
			if genodic[snp] == 'NA':
				continue
			print >> currbimbam, ", ".join(genodic[snp])
		except KeyError:
			#tabixer = pysam.Tabixfile('/mnt/lustre/home/cusanovich/500HT/Imputed1415/ByChr/hutt.imputed.' + chrm + '.txt.gz')
			tabixer = pysam.Tabixfile('/mnt/lustre/home/cusanovich/500HT/' + mapper + '/ByChr/hutt.' + mapper + '.' + distance + '.' + chrm + '.txt.gz')
			for record in tabixer.fetch(chrm,int(snpdic[snp][1]),int(snpdic[snp][2])):
				genos = record.split('\t')
			tabixer.close()
			y = [genos[4], 'A', 'G'] + genos[5:len(genos)]
			missing = len([k for k, j in enumerate(y) if j == 'NA'])
			maf = 1 - (float(len([k for k, j in enumerate(y) if j == '2'])*2 + len([k for k, j in enumerate(y) if j == '1']))/float(len(y)*2 - missing))
			#print maf
			if missing > 21:
				genodic[snp] = 'NA'
				continue
			if maf < 0.05:
				genodic[snp] = 'NA'
				continue
			genodic[snp] = y
			print >> currbimbam, ", ".join(genodic[snp])
	currbimbam.close()
	phener = ('cut -f' + str(int(exprcoldic[gene]) + 1) + ' -d" " ' + hmdir + '500HT/qqnorm.500ht.' + mapper + '_order.bimbam > ' + currfiles + '.pheno')
	ifier(phener)
	#print "Running GEMMA..."
	gemmer = ('cd ' + genodir + '; ' + hmdir + 'Programs/gemma0.94 -g ' + currfiles + '.bimbam -p ' + currfiles + '.pheno -k ' + currfiles + '.square.txt -c ' + currfiles + '.pcs.txt -lmm 2 -o curr_' + chrm + '_pc' + str(pcs))
	ifier(gemmer)
	currresults = open(genodir + '/output/curr_' + chrm + '_pc' + str(pcs) + '.assoc.txt','r')
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
	for perm in range(10000):
		if broken == 1:
			continue
		currperm = open(genodir + 'perm_curr_' + chrm + '_pc' + str(pcs) + '.bimbam','w')
		randind = random.shuffle(range(431))
		for snp in masterdic[gene]:
			print >> currperm, ", ".join(genodic[snp][randind])
		currperm.close()
		print str(winnerperms)
		gemmer = ('cd ' + genodir + '; ' + hmdir + 'Programs/gemma0.94 -g ' + genodir + 'perm_curr_' + chrm + '_pc' + str(pcs) + '.bimbam -p ' + currfiles + '.pheno -k ' + currfiles + '.square.txt -c ' + currfiles + '.pcs.txt' + ' -lmm 2 -o perm_curr_' + chrm + '_pc' + str(pcs))
		ifier(gemmer)
		perms = []
		permer = open(genodir + 'output/perm_curr_' + chrm + '_pc' + str(pcs) + '.assoc.txt','r')
		permlow = 1
		for line in permer:
			liner = line.strip().split()
			if 'nan' in liner[5] or 'p_lrt' in liner[5]:
				continue
			if float(liner[5]) < permlow:
				permlow = float(liner[5])
		permer.close()
		if permlow <= pmin:
			winnerperms += 1
		if winnerperms == 10:
			genep = 11/random.uniform(perm+2,perm+3)
			broken = 1
	if broken == 0:
		genep = (winnerperms + 1)/float(10001)
	winnerdic[gene] = [snpmin, pmin, genep]
	cleanup = 'rm ' + hmdir + '500HT/3chip/*curr_' + chrm + '_pc' + str(pcs) + '*'
	ifier(cleanup)
#	completedgenes += 1

#t1 = time.time()
#print t1-t0

print "Writing results..."
aller = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.3chip.gemma.eqtls.txt','w')
for i in range(len(genes)):
	print >> aller, '{0}\t{1}\t{2:.4g}'.format(genes[i],snps[i],pvals[i])

aller.close()

winners = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.3chip.gemma.chosen.txt','w')
for gene in sorted(winnerdic.keys()):
	print >> winners, '{0}\t{1[0]}\t{1[1]:.4g}\t{1[2]:.4g}'.format(gene,winnerdic[gene])

winners.close()

doners = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.done','w')
doners.close()
