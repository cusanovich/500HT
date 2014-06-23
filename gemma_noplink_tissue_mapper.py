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
from DarrenTools import ifier, matrix_reader
#import time

chrm = sys.argv[1]
pcs = sys.argv[2]
#chrm = 'chr22'
#pcs = 1
genodir = '/mnt/lustre/home/cusanovich/500HT/Imputed1415/'
os.chdir(genodir)
hmdir = '/mnt/lustre/home/cusanovich/'

currfiles = genodir + "curr_" + chrm + "_pc" + str(pcs) + "_" + correction

# This is the bit about running permutations
def permer(scores=currscores, circuit=currperms, permwins=currpermwins, actives=curractive):
	for perm in xrange(0,10000):
		if actives.count(0) == 0:
			continue
		# First we permute the genotypes and write them out to a file (bimbam format)
		shuffle(randind)
		updateind = [0,1,2] + randind
		permbimbam = open(genodir + 'perm_curr.bimbam','w')
		for snp in masterdic[gene]:
			tabixer = Tabixfile('/mnt/lustre/home/cusanovich/500HT/Imputed1415/ByChr/hutt' + mapper + '.' + chrm + '.txt.gz')
			genos = [x.split('\t') for x in tabixer.fetch(chrm,int(snpdic[snp][1]),int(snpdic[snp][2]))][0]
			tabixer.close()
			y = [genos[3], 'A', 'G'] + genos[6:len(genos)]
			yrand = y[updateind]
			print >> permbimbam, ", ".join(yrand)
		permbimbam.close()
		# This runs gemma on the permuted genotypes
		gemmer = (hmdir + 'Programs/gemma0.94 -g ' + genodir + 'perm_curr.bimbam -p ' + currfiles + '.pheno -k ' + currfiles + '.square.txt -c ' + currfiles + '.covs.txt' + ' -lmm 4 -maf 0.05 -o perm_curr')
		ifier(gemmer)
		permresults = matrix_reader(genodir + 'output/perm_curr.assoc.txt',dtype='f8')
		permsort = permresults[permresults[:,12].argsort()]
		permwins = permsort[0:100,]
		permscores = [0]*len(dhsdic[dhs.keys()[0]])
		for snp in permwins[:,1]:
			permscores = permscores + dhsdic[snp]
		# For eQTLs I was picking the smallest cis p-val, we should change this to run the enrichment scripts and then count winners for each tissue
		for z in xrange(len(permscores)):
			if actives[z] == 1:
				continue
			if permscores[z] > 	
		if permlow <= pmin:
			winnerperms += 1
		if winnerperms == 10:
			return 11/uniform(perm+2,perm+3)
	return (winnerperms + 1)/float(10001)

print "Loading expression..."
# This built up a list cis eQTLs, we can remove references to it

overlaps = open(['SOME FILE OF THE SNPs OLAPPING DNASE'],'r')
dnasetable = DarrenTools.matrix_reader(['SAME FILE'])
masterdic = {}
for line in overlaps:
	liner = line.strip().split()
	masterdic[liner[3]] = mastercols[0:3]
	## Also make a dhsdic dictionary of SNPs and the HSs they overlap

#Just copies covariance matrix so parallel processes don't all have to hit the same file over and over
cover = ('cp ' + hmdir + '500HT/addSNP.500ht.ordered.square.txt ' + currfiles + '.square.txt')
ifier(cover)

#Copy phenotype fill to curr directory
cover = ('cp ' + hmdir + '[SOME GWAS PHENOTYPE FILE] ' + currfiles + '.pheno')
ifier(cover)

#completedgenes = 0
pvals = []
genes = []
snps = []
winnerdic = {}
genodic = {}
randind = range(3,434)
#t0 = time.time()

for gene in masterdic.keys():
#	if completedgenes == 5:
#		break
	####Pull genotypes for the SNPs in cis, if genotypes not already in dictionary: go to geno file and pull in appropriate data
	currbimbam = open(currfiles + '.bimbam','w')
	for snp in masterdic[gene]:
		tabixer = Tabixfile('/mnt/lustre/home/cusanovich/500HT/Imputed1415/ByChr/hutt' + mapper + '.' + chrm + '.txt.gz')
		genos = [x.split('\t') for x in tabixer.fetch(chrm,int(snpdic[snp][1]),int(snpdic[snp][2]))][0]
		tabixer.close()
		y = [genos[3], 'A', 'G'] + genos[6:len(genos)]
		print >> currbimbam, ", ".join(y)
	currbimbam.close()
	#print "Running GEMMA..."
	gemmer = (hmdir + 'Programs/gemma0.94 -g ' + currfiles + '.bimbam -p ' + currfiles + '.pheno -k ' + currfiles + '.square.txt -lmm 4 -maf 0.05 -o curr_' + chrm + '_pc' + str(pcs) + '_' + correction)
	ifier(gemmer)
	currresults = open(genodir + 'output/curr_' + chrm + '_pc' + str(pcs) + '_' + correction + '.assoc.txt','r')
	currresults = matrix_reader(genodir + 'output/curr_' + chrm + '_pc' + str(pcs) + '_' + correction + '.assoc.txt',dtype='f8')
	currsort = curresults[currresults[:,12].argsort()]
	currwins = currsort[0:100,]
	currscores = [0]*len(dhsdic[dhs.keys()[0]])
	for snp in currwins[:,1]:
		currscores = currscores + dhsdic[snp]
	currperms = [0]*len(dhsdic[dhs.keys()[0]])
	currpermwins = [0]*len(dhsdic[dhs.keys()[0]])
	curractive = [0]*len(dhsdic[dhs.keys()[0]])
	tissueps = permer(scores = currscores, circuit = currperms, permwins = currpermwins, actives=curractive)
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
