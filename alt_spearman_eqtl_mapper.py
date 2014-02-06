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
import scipy
from scipy import stats
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
mastercols = matrix_reader('/mnt/lustre/home/cusanovich/500HT/hutt.3chip.500kb.mastercols.txt',dtype='|S15')
masterdic = {}
exprcoldic = {}
for i in range(mastercols.shape[0]):
	try:
		masterdic[mastercols[i,0]].append(mastercols[i,1])
	except KeyError:
		masterdic[mastercols[i,0]] = [mastercols[i,1]]
		exprcoldic[mastercols[i,0]] = mastercols[i,2]

naming = open('/mnt/lustre/home/cusanovich/500HT/findivs.500ht.txt','r')
exprnames = naming.readlines()
exprnames = [x.strip().split()[0] for x in exprnames]
exprgenes = list(set(mastercols[:,0]))
exprdic = {}
for gene in exprcoldic.keys():
	exprdic[gene] = exprs[:,exprcoldic[gene]]

print "Loading SNP annotations..."
snpdic = {}
snpbed = open('/mnt/lustre/home/cusanovich/500HT/hutt.3chip.hg19.bed','r')
for line in snpbed:
	liner = line.strip().split()
	snpdic[liner[3]] = liner[0:3]

#genonames = matrix_reader('/mnt/lustre/home/cusanovich/500HT/Imputed1415/imputed_cgi.fam',sep=" ")
genonames = matrix_reader('/mnt/lustre/home/cusanovich/oldhome_nb/Hutterite_Heritability/gemma/Genotypes/hutt.3chip.hg19remap.fam',sep=" ")
genonames = list(genonames[:,1])
genonames = ['chrm','start','end','score','strand'] + genonames
try:
	genoindex = [genonames.index(z) for z in exprnames]
except ValueError:
	print z
	sys.exit()

#1 - for each gene, collect all SNPs
#2 - for each SNP calculate p-value in R
#3 - pick minimum p-value (pick closest if tied)
#4 - permute genotypes (keep LD structure) x 10000
#5 - recalculate p-values and pick mins in R
#6 - after each permutation, check if 10 permutations more extreme -> estimate perm p-value in R
#7 - calculate p-value -> (1 + sum of perms more extreme)/(1 + No. Permutations)
#8 - write out all p-values
#9 - write out gene-level p-values and chosen SNP
obsp = robjects.r('function(exprs,genos){\n'
				 'x = as.numeric(unlist(exprs))\n'
				 'y = as.numeric(unlist(genos))\n'
				 'corr = cor.test(x,y,method="spearman",exact=FALSE)$estimate\n'
				 'pval = cor.test(x,y,method="spearman",exact=FALSE)$p.value\n'
				 'return(c(corr,pval))}')
corp = robjects.r('function(exprs,genos,obsp){\n'
				  'x = as.numeric(unlist(exprs))\n'
				  'y = matrix(unlist(genos),length(unlist(genos))/431,431)\n'
				  'z = as.numeric(unlist(obsp))\n'
				  'winners = 0\n'
				  'for(i in 1:10000){\n'
				  'currgenos = y[,sample(431)]\n'
				  'if(dim(y)[1] == 1){permp = cor.test(as.numeric(y),x,method="spearman",exact=FALSE)$p.value\n'
				  '}else{permp = min(apply(currgenos,1,function(b){cor.test(as.numeric(b),x,method="spearman",exact=FALSE)$p.value}))}\n'
				  'if(permp =< z){winners = winners + 1}\n'
				  'if(winners == 10){\n'
				  #'print(paste0("Killed loop after ",i," permutations."))\n'
				  'return(11/runif(1,i+1,i+2))}}\n'
				  'return((winners + 1)/10001)}')

obs = []
pvals = []
genes = []
snps = []
minors = []
genodic = {}
winnerdic = {}
#t0 = time.time()
for gene in masterdic.keys():
#for gene in masterdic.keys()[1:500]:
	if snpdic[masterdic[gene][0]][0] != chrm:
		continue
	x = exprdic[gene]
	snping = [0]*len(masterdic[gene])
	snpmin = []
	genemin = []
	pmin = 1
	genolist = []
	for i, snp in enumerate(masterdic[gene]):
		try:
			y = genodic[snp]
			if genodic[snp] == 'NA':
				continue
		except KeyError:
			#tabixer = pysam.Tabixfile('/mnt/lustre/home/cusanovich/500HT/Imputed1415/ByChr/hutt.imputed.' + chrm + '.txt.gz')
			tabixer = pysam.Tabixfile('/mnt/lustre/home/cusanovich/oldhome_nb/Hutterite_Heritability/gemma/Genotypes/ByChr/hutt.3chip.' + chrm + '.txt.gz')
			for record in tabixer.fetch(chrm,int(snpdic[snp][1]),int(snpdic[snp][2])):
				genos = record.split('\t')
			tabixer.close()
			y = [genos[index] for index in genoindex]
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
		currcor = obsp(list(x),y)
		snping[i] = currcor[1]
		#print currcor[1]
		obs.append(currcor[0])
		pvals.append(currcor[1])
		genes.append(gene)
		snps.append(snp)
		genolist = genolist + y
		if currcor[1] < pmin:
			snpmin = snp
			pmin = currcor[1]
	genep = corp(list(x),genolist,pmin)
	winnerdic[gene] = [snpmin, pmin, genep[0]]

#t1 = time.time()
#print t1-t0

print "Writing results..."
aller = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.3chip.500kb.eqtls.txt','w')
for i in range(len(genes)):
	print >> aller, '{0}\t{1}\t{2:.3g}\t{3:.4g}'.format(genes[i],snps[i],obs[i],pvals[i])

aller.close()

winners = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.3chip.500kb.chosen.txt','w')
for gene in sorted(winnerdic.keys()):
	print >> winners, '{0}\t{1[0]}\t{1[1]:.4g}\t{1[2]:.4g}'.format(gene,winnerdic[gene])

winners.close()

doners = open('/mnt/lustre/home/cusanovich/500HT/ByChr/' + chrm + '.PC' + str(pcs) + '.done','w')
doners.close()