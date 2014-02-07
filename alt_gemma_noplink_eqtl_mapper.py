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
time.sleep(random.randint(1,150))
mastercols = matrix_reader('/mnt/lustre/home/cusanovich/500HT/hutt.3chip.mastercols.txt',dtype='|S15')
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
			tabixer = pysam.Tabixfile('/mnt/lustre/home/cusanovich/500HT/3chip/ByChr/hutt.3chip.' + chrm + '.txt.gz')
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
	phener = ('cut -f' + str(int(exprcoldic[gene]) + 1) + ' -d" " ' + hmdir +
		'500HT/qqnorm.500ht.3chip_order.bimbam > ' + currfiles + '.pheno; paste -d" " ' + currfiles + '.header.fam ' + currfiles + '.tailer.fam > ' +
		currfiles + '.fam; rm ' + currfiles + '.header.fam; rm ' + currfiles + '.tailer.fam')
	#print "Running GEMMA..."
	gemmer = ('cd ' + genodir + '; ' + hmdir + 'Programs/gemma0.94 -g ' + currfiles + '.bimbam -p ' + currfiles + '.pheno -k ' + hmdir + '500HT/addSNP.500ht.3chip_order.square.txt -c ' + hmdir + '500HT/Exprs/qqnorm.500ht.3chip_order.pc' + str(pcs) + ' -lmm 2 -n ' + str(int(exprcoldic[gene]) + 1) + ' -o curr_' + chrm + '_pc' + str(pcs))
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
	findic = {}
	for perm in range(10000):
		if broken == 1:
			continue
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
		random.shuffle(randfindivs)
		rearrange = open(currfiles + '.findivs.rearrange.txt','w')
		for z in range(len(findivs)):
			print >> rearrange, 'HUTTERITES\t' + findivs[z] + '\tHUTTERITES\t' + randfindivs[z]
		rearrange.close()
		plinker = 'plink --noweb --bfile ' + currfiles + ' --update-ids ' + currfiles + '.findivs.rearrange.txt --make-bed --out ' + genodir + 'perm_curr_' + chrm + '_pc' + str(pcs)
		ifier(plinker)
		permfam = open(genodir + 'perm_curr_' + chrm + '_pc' + str(pcs) + '.fam','w')
		permpc = open(genodir + 'perm_curr_' + chrm + '_pc' + str(pcs) + '.pcs','w')
		for z in range(len(findivs)):
			print >> permfam, " ".join(findic[randfindivs[z]])
			print >> permpc, " ".join(pcdic[randfindivs[z]])
		permfam.close()
		permpc.close()
		print str(winnerperms)
		gemmer = ('cd ' + genodir + '; ' + hmdir + 'Programs/gemma0.94 -bfile perm_curr_' + chrm + '_pc' + str(pcs) + ' -km 2 -k ' + hmdir + '500HT/addSNP.500ht.txt -c ' + genodir + 'perm_curr_' + chrm + '_pc' + str(pcs) + '.pcs' + ' -lmm 2 -o perm_curr_' + chrm + '_pc' + str(pcs))
		ifier(gemmer)
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
