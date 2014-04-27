#!/usr/bin/env python
import os
import sys
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
import pysam
from pysam import Tabixfile
import subprocess
import numpy
import random
from random import shuffle
from random import uniform
import gzip
#import time

chrm = sys.argv[1]
pcs = sys.argv[2]
#chrm = 'chr22'
#pcs = 1
genodir = '/mnt/lustre/home/cusanovich/500HT/Imputed1415/'
os.chdir(genodir)
hmdir = '/mnt/lustre/home/cusanovich/'

currfiles = genodir + "curr_" + chrm + "_pc" + str(pcs) + "_" + correction

# This function just makes it easier to communicate with the shell from Python
def ifier(commander):
	ify = subprocess.Popen(commander,shell=True)
	ify.wait()

# This function just makes python read in big data files all at once more efficiently
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

# This is the bit about running permutations
def permer(gene):
	winnerperms = 0
	for perm in xrange(0,10000):
		# First we permute the genotypes and write them out to a file (bimbam format)
		shuffle(randind)
		updateind = [0,1,2] + randind
		permgenos = [", ".join([genodic[x][index] for index in updateind]) for x in masterdic[gene]]	#this permutes genotypes
		currbimbam = open(genodir + 'perm_curr_' + chrm + '_pc' + str(pcs) + '.bimbam','w')
		print >> currbimbam, "\n".join(permgenos)
		currbimbam.close()
		# These just report where we are.  We can switch all the gene count to phenotype count.
		print 'Gene No. ' + str(len(winnerdic.keys()) + 1) + ' on chrm.'
		print str(winnerperms) + ' of ' + str(perm) + ' permutations lost.'
		# This runs gemma on the permuted genotypes
		gemmer = (hmdir + 'Programs/gemma0.94 -g ' + genodir + 'perm_curr_' + chrm + '_pc' + str(pcs) + '.bimbam -p ' + currfiles + '.pheno -k ' + currfiles + '.square.txt -c ' + currfiles + '.pcs.txt' + ' -lmm 4 -maf 0.05 -o perm_curr_' + chrm + '_pc' + str(pcs))
		ifier(gemmer)
		# This reads in the GEMMA results.
		permering = open(genodir + 'output/perm_curr_' + chrm + '_pc' + str(pcs) + '.assoc.txt','r')
		permers = [x.strip().split()[12] for x in permering.readlines()]
		permers = [float(x) for x in permers if x != 'nan' and x != 'p_lrt']
		permering.close()
		# For eQTLs I was picking the smallest cis p-val, we should change this to run the enrichment scripts and then count winners for each tissue
		permlow = min(permers)
		if permlow <= pmin:
			winnerperms += 1
		if winnerperms == 10:
			return 11/uniform(perm+2,perm+3)
	return (winnerperms + 1)/float(10001)

print "Loading expression..."
# This built up a list cis eQTLs, we can remove references to it
mastercols = matrix_reader('/mnt/lustre/home/cusanovich/500HT/hutt' + mapper + '.' + distance + '.mastercols.txt',dtype='|S15')

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
snpbed = open('/mnt/lustre/home/cusanovich/500HT/hutt.imputed.coord.bed','r')
for line in snpbed:
	liner = line.strip().split()
	snpdic[liner[3]] = liner[0:3]

#Just copies covariance matrix so parallel processes don't all have to hit the same file over and over
cover = ('cp ' + hmdir + '500HT/addSNP.500ht.ordered.square.txt ' + currfiles + '.square.txt')
ifier(cover)

#Copy phenotype fill to curr directory
cover = ('cp ' + hmdir + '[SOME GWAS PHENOTYPE FILE] ' + currfiles + '.pheno')
ifier(cover)

#chr22ers = []
#for gene in masterdic.keys():
#	if chrmdic[gene] == chrm:
#		chr22ers.append(gene)

#completedgenes = 0
pvals = []
genes = []
snps = []
winnerdic = {}
genodic = {}
randind = range(3,434)
#Time-savers
geneappend = genes.append
snpappend = snps.append
pvalappend = pvals.append
#t0 = time.time()

for gene in masterdic.keys():
	if chrmdic[gene] != chrm:
		continue
#	if completedgenes == 5:
#		break
	print gene
	currgenos = []
	####Pull genotypes for the SNPs in cis, if genotypes not already in dictionary: go to geno file and pull in appropriate data
	for snp in masterdic[gene]:
		try:
			currgenos.append(", ".join(genodic[snp]))
		except KeyError:
			#tabixer = pysam.Tabixfile('/mnt/lustre/home/cusanovich/500HT/Imputed1415/ByChr/hutt.imputed.' + chrm + '.txt.gz')
			#tabixer = pysam.Tabixfile('/mnt/lustre/home/cusanovich/500HT/' + mapper + '/ByChr/hutt.' + mapper + '.' + distance + '.' + chrm + '.txt.gz')
			tabixer = Tabixfile('/mnt/lustre/home/cusanovich/500HT/Imputed1415/ByChr/hutt' + mapper + '.' + chrm + '.txt.gz')
			genos = [x.split('\t') for x in tabixer.fetch(chrm,int(snpdic[snp][1]),int(snpdic[snp][2]))][0]
			tabixer.close()
			y = [genos[3], 'A', 'G'] + genos[6:len(genos)]
			genodic[snp] = y
			currgenos.append(", ".join(genodic[snp]))
	currbimbam = open(currfiles + '.bimbam','w')
	print >> currbimbam, "\n".join(currgenos)
	currbimbam.close()
	#print "Running GEMMA..."
	if regressPCs:
		gemmer = (hmdir + 'Programs/gemma0.94 -g ' + currfiles + '.bimbam -p ' + currfiles + '.pheno -k ' + currfiles + '.square.txt -lmm 4 -maf 0.05 -o curr_' + chrm + '_pc' + str(pcs) + '_' + correction)
		ifier(gemmer)
	if not regressPCs:
		gemmer = (hmdir + 'Programs/gemma0.94 -g ' + currfiles + '.bimbam -p ' + currfiles + '.pheno -k ' + currfiles + '.square.txt -c ' + currfiles + '.pcs.txt -lmm 4 -maf 0.05 -o curr_' + chrm + '_pc' + str(pcs) + '_' + correction)
		ifier(gemmer)
	currresults = open(genodir + '/output/curr_' + chrm + '_pc' + str(pcs) + '_' + correction + '.assoc.txt','r')
	pmin = 1.1
	snpmin = ''
	currcount = 0
	for line in currresults:
		if 'p_lrt' in line or 'nan' in line:
			continue
		liner = line.strip().split()
		geneappend(gene)
		snpappend(liner[1])
		pvalappend(float(liner[12]))
		currcount = currcount + 1
		if float(liner[12]) < pmin:
			pmin = float(liner[12])
			snpmin = liner[1]
	if pmin == 1.1:
		continue
	#print "Starting up permutations..."
	if bonferroni:
		genep = min(1,pmin*currcount)
	else:
		genep = permer(gene)
	winnerdic[gene] = [snpmin, pmin, genep]
	#completedgenes += 1

#t1 = time.time()
#print t1-t0

cleanup = 'rm ' + genodir + '*curr_' + chrm + '_pc' + str(pcs) + '_' + correction + '.*'
ifier(cleanup)

#print "Writing results..."
aller = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + mapper + '.' + distance + '.' + correction + gccor + covcor + regressed + '.gemma.eqtls.txt','w')
for i in xrange(0,len(genes)):
	print >> aller, '{0}\t{1}\t{2:.4g}'.format(genes[i],snps[i],pvals[i])

aller.close()

winners = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + mapper + '.' + distance + '.' + correction + gccor + covcor + regressed + '.gemma.chosen.txt','w')
for gene in sorted(winnerdic.keys()):
	print >> winners, '{0}\t{1[0]}\t{1[1]:.4g}\t{1[2]:.4g}'.format(gene,winnerdic[gene])

winners.close()

doners = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.' + correction + '.done','w')
doners.close()
