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
from DarrenTools import ifier, matrix_reader
#import time

phenofile = 'absneutroCombSNP.ph-cvt'
blocknum = '00'
#phenofile = sys.argv[1]
#blocknum = sys.argv[2]
pheno = phenofile.split('.')[0]
genodir = '/mnt/lustre/home/cusanovich/500HT/Imputed1415/'
os.chdir(genodir + 'perms/')
hmdir = '/mnt/lustre/home/cusanovich/'

currfiles = genodir + "perms/Block_" + blocknum + "_curr_" + pheno

# This is the bit about running permutations
def permer(perm,randind):
		# First we permute the genotypes and write them out to a file (bimbam format)
		random.seed(perm + (int(blocknum)*100))
		shuffle(randind)
		snper = currfiles + '.snps'
		aer = currfiles + '.as'
		ger = currfiles + '.gs'
		shuffler = 'zcat ' + genodir + 'ByChr/*.all*.gz | grep -f ' + snper + ' - | cut -f' + ','.join(randind) + ' > ' + currfiles + '_perm_sub.bimbam; paste ' + snper + ' ' + aer + ' ' + ger + ' ' + currfiles + '_perm_sub.bimbam > ' + currfiles + '_perm.bimbam'
		ifier(shuffler)
		# This runs gemma on the permuted genotypes
		gemmer = (hmdir + 'Programs/gemma0.94 -g ' + currfiles + '_perm.bimbam -p ' + currfiles + '.pheno -k ' + currfiles + '.square.txt -c ' + currfiles + '.covariates.txt' + ' -lmm 4 -maf 0.05 -o ' + blocknum + '.' + pheno)
		ifier(gemmer)

print "Loading data..."
overlaps = open(genodir + 'hutt.imputed.dnase2.bed','r')
#dnasetable = DarrenTools.matrix_reader(genodir + 'hutt.imputed.dnase.bed')
masterdic = {}
for line in overlaps:
	liner = line.strip().split()
	if liner[4:].count('0') == len(liner[4:]):
		continue
	masterdic[liner[3]] = liner[0:3]

overlaps.close()
exclfile = open(hmdir + '500HT/qqnorm.500ht.gccor.newcovcor.findivs.txt','r')
exclusions = exclfile.readlines()
exclusions = [x.strip('\n') for x in exclusions]
exclfile.close()
#Just copies covariance matrix so parallel processes don't all have to hit the same file over and over

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

####Pull genotypes for the SNPs in cis, if genotypes not already in dictionary: go to geno file and pull in appropriate data
print "Writing out genotypes..."
genofinfile = open(genodir + 'hutt.imputed.rename.fam','r')
genofins = []
for line in genofinfile:
	genofins.append(line.strip().split()[1])

genofinfile.close()

genoinds = [genofins.index(x) + 6 for x in officialfindivs]

snplist = open(currfiles + '.snps','w')
alist = open(currfiles + '.as','w')
glist = open(currfiles + '.gs','w')
for snp in masterdic.keys():
#for snp in masterdic.keys()[0:1000]:
	print >> snplist, snp
	print >> alist, 'A'
	print >> glist, 'G'

snplist.close()
alist.close()
glist.close()

print "Running permutations..."
blocker = open(genodir + 'Block_' + blocknum + 'permwins.txt','w')
for perm in xrange(0,100):
	permer(perm,randind = genoinds)
	permresults = matrix_reader(genodir + 'output/perm_curr_' + blocknum + '.assoc.txt',dtype='f8')
	permsort = permresults[permresults[:,12].argsort()]
	permwins = permsort[0:100,]
	print >> blocker, '\t'.join(permwins)

blocker.close()
 
cleanup = 'rm ' + genodir + '*curr_' + pheno + '.*'
ifier(cleanup)

print "Writing results..."
aller = open('/mnt/lustre/home/cusanovich/500HT/Tissues/' + pheno + '.enrichmentps.txt','w')
for i in xrange(0,len(tissueps)):
	print >> aller, '{0}\t{2:.4g}'.format(dhsdic['rsID'][i],tissueps[i])

aller.close()
