#!/usr/bin/env python
import os
import sys
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
from myfuncs import matrix_reader,ifier
import subprocess
import gzip
import numpy as np
#import time
#import pysam
import rpy2.robjects as robjects
#import random

#1 - iterate through pval files
#2 - For each row:
#		count number of cis occurences
#		choose smallest p and bonferroni correct and record snp and gene names
#3 - write file with gene name snp name original p corrected p
#4 - iterate back through all files and re-collect all ps that match correct global bonferroni and write those out to a file

print "Loading annotations..."
sys.stdout.flush()
qtldir = '/mnt/lustre/home/cusanovich/500HT/ByChr/'
os.chdir(qtldir)
files = os.listdir(qtldir)
files = [x for x in files if 'trans.pvals' in x]
#chrm = '22'
#Number of genes * Number of SNPs less any cis associations
#bonthresh = 14014*393613.0

rsfile = open('/mnt/lustre/home/cusanovich/500HT/Imputed1415/trans.rsids.txt','r')
rsnos = [x.strip().split()[0] for x in rsfile.readlines()]
rsfile.close()

tssdic = {}
tss = open('../ensemblCAGETSS_RNAseq_sorted.bed','r')
for line in tss:
	liner = line.strip().split()
	tssdic[liner[3]] = liner[1]

snps = open('../Imputed1415/hutt.imputed.500ht.bim','r')
snpdic = {}
for line in snps:
	liner = line.strip().split()
	snpdic[liner[1]] = ['chr' + liner[0], liner[3]]

print "Collecting cislike snps..."
sys.stdout.flush()
notcisers = 0
winners = []
winningps = []
winningbonf = []
genes = []
sherlockfile = open('../eQTLs/remaster.imputed.1Mb.bonferroni.gccor.newcovcor.sherlock.txt','w')
for chrm in files:
	currchrm = gzip.open(chrm)
	for line in currchrm:
		#if len(genes) >= 50:
		#	continue
		currnotcisers = 0
		currwinner = ''
		currwinning = 1
		liner = line.strip().split()
		genes.append(liner[0])
		print "Gene no. " + str(len(genes)) + " (" + chrm + ")"
		sys.stdout.flush()
		for counter,snp in enumerate(liner):
			try:
				tester = float(snp)
				currnotcisers += 1
				if tester < 0.001:
					print >> sherlockfile, liner[0] + '\t' + rsnos[counter] + '\t' + str(tester) + '\t0'
				if tester < currwinning:
					currwinning = float(snp)
					currwinner = rsnos[counter]
			except ValueError:
				if 'ENSG' in snp:
					continue
				try:
					tester = float(snp.split('(cis)')[0])
					if tester < 0.1:
						print >> sherlockfile, liner[0] + '\t' + rsnos[counter] + '\t' + str(tester) + '\t1'
				except ValueError:
					print 'No. ' + str(counter) + ': ' + chrm + ' ' + rsnos[counter] + ' bad!'
					pass
		winners.append(currwinner)
		winningps.append(currwinning)
		winningbonf.append(min(currwinning*currnotcisers,1))
		notcisers = notcisers + currnotcisers

sherlockfile.close()
print "Writing cislike file..."
sys.stdout.flush()
cislike = open('../eQTLs/master.imputed.1Mb.bonferroni.gccor.newcovcor.trans.cislike.pvals.txt','w')
for counter,gene in enumerate(genes):
	newestline = [gene, winners[counter], str(winningps[counter]), str(winningbonf[counter])]
	print >> cislike, '\t'.join(newestline)

cislike.close()

bonthresh = float(notcisers)
print "Bonferroni threshold is " + str(1/bonthresh)
print "Writing trans file..."
genecount = 0
sys.stdout.flush()
transers = open('../eQTLs/master.imputed.1Mb.bonferroni.gccor.newcovcor.trans.sig.pvals.txt','w')
for chrmagain in files:
	currchrm = gzip.open(chrmagain)
	chrmed = chrmagain.split('.')[0]
	for line in currchrm:
		#if genecount >= 50:
		#	continue
		genecount += 1
		liner = line.strip().split()
		print "Gene no. " + str(genecount) + " again (" + chrm + ")"
		sys.stdout.flush()
		for counter,snp in enumerate(liner):
			try:
				tester = float(snp)
				if tester < 1/bonthresh:
					snpinfo = snpdic[rsnos[counter]]
					if snpinfo[0] == chrmed and abs(int(snpinfo[1]) - int(tssdic[liner[0]])) < 1000000:
						continue
					print >> transers, rsnos[(counter + 1)] + '\t' + snpinfo[0] + ':' + snpinfo[1] + '\t' + liner[0] + '\t' + chrmed + ':' + tssdic[liner[0]] + '\t' + snp
			except ValueError:
				pass

transers.close()

