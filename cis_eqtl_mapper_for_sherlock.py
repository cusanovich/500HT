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

#if len(sys.argv) != 3:
#	print "Usage = python gemma_eqtl_mapper.py chr pc"
#	sys.exit()
gccorrection = True
covcorrection = True
bonferroni = True
regressPCs = True
chrm = sys.argv[1]
pcs = sys.argv[2]
#chrm = 'chr22'
#pcs = 1
genodir = '/mnt/lustre/home/cusanovich/500HT/Imputed1415/'
os.chdir(genodir)
hmdir = '/mnt/lustre/home/cusanovich/'
mapper = '.imputed'
distance = '150kb'
if bonferroni:
	correction = 'bonferroni'
else:
	correction = 'permutation'

gccor = ''
covcor = ''
regressed = ''
if gccorrection:
	gccor = '.gccor'

if covcorrection:
	covcor = '.covcor'

if regressPCs:
	regressed = '.regressPCs'

currfiles = genodir + "curr_" + chrm + "_pc" + str(pcs) + "_" + correction

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
if regressPCs:
	mastercols = matrix_reader('/mnt/lustre/home/cusanovich/500HT/hutt' + mapper + '.' + distance + '.mastercols.txt',dtype='|S15')

if not regressPCs:
	mastercols = matrix_reader('/mnt/lustre/home/cusanovich/500HT/hutt' + mapper + '.' + distance + '.chrmspecific.mastercols.txt',dtype='|S15')

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

cover = ('cp ' + hmdir + '500HT/addSNP.500ht.ordered.square.txt ' + currfiles + '.square.txt')
ifier(cover)

if regressPCs:
	expressed = matrix_reader(hmdir + '500HT/qqnorm.500ht' + gccor + covcor + '.ordered.bimbam',dtype='float')
	pcmat = matrix_reader(hmdir + '500HT/Exprs/qqnorm.500ht' + gccor + covcor + '.ordered.pc' + str(pcs),dtype='float')

if not regressPCs:
	cover = ('cp ' + hmdir + '500HT/Exprs/qqnorm.500ht' + gccor + covcor + '.ordered.pc' + str(pcs) + ' ' + currfiles + '.pcs.txt')
	ifier(cover)

if regressPCs and int(pcs) != 0:
	mod1 = pcmat
	mod2 = mod1.T
	Y = expressed.T
	W = pcmat[:,1:(int(pcs)+1)]
	mods = mod2.dot(mod1)
	inversemods = numpy.linalg.inv(mods)
	gammahat = Y.dot(mod1).dot(inversemods)[:,1:mods.shape[1]]
	Yfit = Y - gammahat.dot(W.T)

if regressPCs and int(pcs) == 0:
	Yfit = expressed.T

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
	if regressPCs:
		numpy.savetxt(currfiles + '.pheno',Yfit[exprcoldic[gene],],delimiter='\n',fmt='%s')
	if not regressPCs:
		phener = ('cut -f' + str(int(exprcoldic[gene]) + 1) + ' -d" " ' +
			hmdir + '500HT/Exprs/qqnorm.500ht' + gccor + covcor +
			'.ordered.' + chrm + '.bimbam > ' + currfiles + '.pheno')
		ifier(phener)
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
