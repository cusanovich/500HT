#!/usr/bin/env python
import os
import sys
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
import subprocess
import numpy
import pysam
import gzip
import rpy2.robjects as robjects
#import time

if len(sys.argv) != 3:
	print "Usage = python spearman_eqtl_mapper.py chr pc"
	sys.exit()

chrm = sys.argv[1]
pcs = sys.argv[2]
#chrm = 'chr22'

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
exprs = matrix_reader('/mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.bimbam.PC' + str(pcs) + '.fixed',dtype='|S10')
mastercols = matrix_reader('/mnt/lustre/home/cusanovich/500HT/hutt.imputed.mastercols.txt',dtype='|S10')
naming = open('/mnt/lustre/home/cusanovich/500HT/findivs.500ht.txt','r')
exprnames = naming.readlines()
exprnames = [x.strip().split()[0] for x in exprnames]
exprgenes = list(set(mastercols[:,0]))
exprdic = {}
for i in range(len(exprgenes)):
	exprdic[exprgenes[i]] = exprs[:,mastercols[i,2]]

print "Loading SNP annotations..."
snpdic = {}
snpbed = open('/mnt/lustre/home/cusanovich/500HT/hutt.imputed.bed','r')
for line in snpbed:
	liner = line.strip().split()
	snpdic[liner[3]] = liner[0:3]

genonames = matrix_reader('/mnt/lustre/home/cusanovich/500HT/Imputed1415/imputed_cgi.fam',sep=" ")
genonames = list(genonames[:,1])
genonames = ['chrm','start','end','score','strand'] + genonames
try:
	genoindex = [genonames.index(z) for z in exprnames]
except ValueError:
	print z
	sys.exit()
corp = robjects.r('function(x,y){\n'
				  'x = as.numeric(unlist(x))\n'
				  'y = as.numeric(unlist(y))\n'
				  'cor0 = cor(x,y,method="spearman",use="pairwise.complete.obs")\n'
				  'z = replicate(1000,cor(x,sample(y),method="spearman",use="pairwise.complete.obs"))\n'
				  'pval = mean(abs(z) >= abs(cor0))\n'
				  'if(pval < 0.01){\n'
				  'zprime = replicate(9000,cor(x,sample(y),method="spearman",use="pairwise.complete.obs"))\n'
				  'pval = mean(abs(c(z,zprime)) >= abs(cor0))}\n'
				  'return(c(cor0,pval))}')

obs = []
pvals = []
genes = []
snps = []
minors = []
genodic = {}
winnerdic = {}
#t0 = time.time()
for i in range(mastercols.shape[0]):
#for i in range(100000,200000):
	if snpdic[mastercols[i,1]][0] != chrm:
		continue
	try:
		if genodic[mastercols[i,1]] == 'NA':
			continue
	except KeyError:
		pass
	x = exprdic[mastercols[i,0]]
	try:
		y = genodic[mastercols[i,1]]
	except KeyError:
		tabixer = pysam.Tabixfile('/mnt/lustre/home/cusanovich/500HT/Imputed1415/ByChr/hutt.imputed.' + chrm + '.txt.gz')
		for record in tabixer.fetch(chrm,int(snpdic[mastercols[i,1]][1]),int(snpdic[mastercols[i,1]][2])):
			genos = record.split('\t')
		tabixer.close()
		y = [genos[index] for index in genoindex]
		missing = len([k for k, j in enumerate(y) if j == 'NA'])
		maf = 1 - (float(len([k for k, j in enumerate(y) if j == '2'])*2 + len([k for k, j in enumerate(y) if j == '1']))/float(860 - missing))
		#print maf
		if missing > 21:
			genodic[mastercols[i,1]] = 'NA'
			continue
		if maf < 0.05:
			genodic[mastercols[i,1]] = 'NA'
			continue
		genodic[mastercols[i,1]] = y
	currcor = corp(list(x),y)
	#print currcor[1]
	obs.append(currcor[0])
	pvals.append(currcor[1])
	genes.append(mastercols[i,0])
	snps.append(mastercols[i,1])
	try:
		winnerdic[mastercols[i,0]][2] = min(currcor[1],winnerdic[mastercols[i,0]][2])
		if currcor[1] == winnerdic[mastercols[i,0]][2]:
			winnerdic[mastercols[i,0]][1] = currcor[0]
			winnerdic[mastercols[i,0]][0] = mastercols[i,1]
	except KeyError:
		winnerdic[mastercols[i,0]] = [mastercols[i,1],currcor[0],currcor[1]]

#t1 = time.time()
#print t1-t0

print "Writing results..."
aller = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.imputed.eqtls.txt','w')
for i in range(len(genes)):
	print >> aller, '{0}\t{1}\t{2:.3f}\t{3:.4f}'.format(genes[i],snps[i],obs[i],pvals[i])

aller.close()

winners = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.imputed.chosen.txt','w')
for gene in sorted(winnerdic.keys()):
	print >> winners, '{0}\t{1[0]}\t{1[1]:.3f}\t{1[2]:.4f}'.format(gene,winnerdic[gene])

winners.close()

doners = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.done','w')
doners.close()