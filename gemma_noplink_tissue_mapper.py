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

phenofile = 'absneutroCombSNP.ph-cvt'
#phenofile = sys.argv[1]
pheno = phenofile.split('.')[0]
genodir = '/mnt/lustre/home/cusanovich/500HT/Imputed1415/'
os.chdir(genodir)
hmdir = '/mnt/lustre/home/cusanovich/'

currfiles = genodir + "curr_" + pheno

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
		permscores = [0]*len(dhsdic[dhsdic.keys()[0]])
		for snp in permwins[:,1]:
			permscores = permscores + dhsdic[snp]
		# For eQTLs I was picking the smallest cis p-val, we should change this to run the enrichment scripts and then count winners for each tissue
		for z in xrange(len(permscores)):
			if actives[z] == 1:
				continue
			if permscores[z] > scores[z]:
				permwins[z] += 1
				circuit[z] += 1
				if permwins[z] == 10:
					actives[z] = 1
	ps = [0]*len(scores)
	for a in xrange(len(permscores)):
		if circuit[a] < 10000:
			ps[a] = 11/uniform(circuit[a]+2,circuit[a]+3)
		else:
			ps[a] = (permwins[a] + 1)/float(circuit[a])
	return ps

print "Loading data..."
overlaps = open(genodir + 'hutt.imputed.dnase.bed','r')
#dnasetable = DarrenTools.matrix_reader(genodir + 'hutt.imputed.dnase.bed')
masterdic = {}
dhsdic = {}
for line in overlaps:
	liner = line.strip().split()
	if liner[4:].count('0') == len(liner[4:]):
		continue
	masterdic[liner[3]] = liner[0:3]
	if liner[3] == 'rsID':
		dhsdic[liner[3]] = liner[4:]
	else:
		dhsdic[liner[3]] = [int(x) for x in liner[4:]]

overlaps.close()
exclfile = open(hmdir + '500HT/qqnorm.500ht.gccor.newcovcor.findivs.txt','r')
exclusions = exclfile.readlines()
exclusions = [x.strip('\n') for x in exclusions]
exclfile.close()
#Just copies covariance matrix so parallel processes don't all have to hit the same file over and over

#cover = ('cp ' + hmdir + '500HT/addSNP.500ht.ordered.square.txt ' + currfiles + '.square.txt')
#ifier(cover)

#Copy phenotype file to curr directory
officialfindivs = []
phener = open(hmdir + '500HT/dege/' + phenofile,'r')
currphen = open(currfiles + '.pheno','w')
currcov = open(currfiles + '.covariates','w')
for line in phener:
	if 'findiv' in line:
		continue
	liner = line.strip().split()
	if liner[0] in exclusions:
		continue
	print >> currphen, liner[1]
	print >> currcov, ' '.join(liner[2:])
	officialfindivs.append(liner[0])

phener.close()
currphen.close()
currcov.close()

covmat = numpy.zeros((len(officialfindivs),len(officialfindivs)))
covs = open(hmdir + 'oldhome_nb/Hutterite_Heritability/Idcoefs/addSNP.coef.3671','r')
for line in covs:
	if 'Ind1' in line:
		continue
	liner = line.strip().split()
	if liner[0] not in officialfindivs or liner[1] not in officialfindivs:
		continue
	ind1 = officialfindivs.index(liner[0])
	ind2 = officialfindivs.index(liner[1])
	covmat[ind1,ind2] = liner[2]
	covmat[ind2,ind1] = liner[2]

covs.close()
numpy.savetxt(currfiles + '.square.txt', covmat, fmt = '%s', delimiter = '\t', newline = '\n')

#Set up list for permutations
randind = range(3,len(officialfindivs) + 3)

#t0 = time.time()
####Pull genotypes for the SNPs in cis, if genotypes not already in dictionary: go to geno file and pull in appropriate data
print "Writing out genotypes..."
genofinfile = open(genodir + 'hutt.imputed.rename.fam','r')
genofins = []
for line in genofinfile:
	genofins.append(line.strip().split()[1])

genofinfile.close()

genoinds = [genofins.index(x) + 6 for x in officialfindivs]
genos = {}
currbimbam = open(currfiles + '.bimbam','w')
#t0 = time.time()
for snp in masterdic.keys():
#for snp in masterdic.keys()[0:1000]:
	chrm = masterdic[snp][0]
	if chrm == 'chrm':
		continue
	tabixer = Tabixfile('/mnt/lustre/home/cusanovich/500HT/Imputed1415/ByChr/hutt.all.imputed.' + chrm + '.txt.gz')
	tempgenos = [x.split('\t') for x in tabixer.fetch(chrm,int(masterdic[snp][1])-1,int(masterdic[snp][2]))][0]
	genos[snp] = [tempgenos[x] for x in range(0,6) + genoinds]
	tabixer.close()
	y = [genos[snp][3], 'A', 'G'] + genos[snp][6:]
	print >> currbimbam, ", ".join(y)

#t1 = time.time()
#print t1-t0
currbimbam.close()

#genomat = matrix_reader(genodir + 'hutt.imputed.dhssnps.bimbam',sep=",")
print "Running GEMMA..."
gemmer = (hmdir + 'Programs/gemma0.94 -g ' + currfiles + '.bimbam -p ' + currfiles + '.pheno -k ' + currfiles + '.square.txt -c ' + currfiles + '.covariates -lmm 4 -maf 0.05 -o curr_' + pheno)
t0 = time.time()
ifier(gemmer)
t1 = time.time()
print t1-t0
#currresults = open(genodir + 'output/curr_' + pheno + '.assoc.txt','r')
currresults = matrix_reader(genodir + 'output/curr_' + pheno + '.assoc.txt',dtype='f8')
currsort = curresults[currresults[:,12].argsort()]
currwins = currsort[0:100,]
currscores = [0]*len(dhsdic[dhsdic.keys()[0]])
for snp in currwins[:,1]:
	currscores = currscores + dhsdic[snp]
currperms = [0]*len(dhsdic[dhsdic.keys()[0]])
currpermwins = [0]*len(dhsdic[dhsdic.keys()[0]])
curractive = [0]*len(dhsdic[dhsdic.keys()[0]])

print "Running permutations..."
tissueps = permer(scores = currscores, circuit = currperms, permwins = currpermwins, actives=curractive)
#completedgenes += 1

#t1 = time.time()
#print t1-t0

cleanup = 'rm ' + genodir + '*curr_' + pheno + '.*'
ifier(cleanup)

print "Writing results..."
aller = open('/mnt/lustre/home/cusanovich/500HT/Tissues/' + pheno + '.enrichmentps.txt','w')
for i in xrange(0,len(tissueps)):
	print >> aller, '{0}\t{2:.4g}'.format(dhsdic['rsID'][i],tissueps[i])

aller.close()
