#!/usr/bin/env python

import sys
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import numpy
import subprocess

biomart = importr("biomaRt")
snpmart = biomart.useMart("snp", dataset="hsapiens_snp")
chrmreq = sys.argv[1]

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

complementers = {'A':'T', 'C':'G', 'T':'A', 'G':'C', '-':'-'}
def reverse_complementer(x):
	x = x.upper()
	splitup = list(x)
	rev = splitup[::-1]
	revcomp = [complementers[z] for z in rev]
	revcomp = ''.join(revcomp)
	return revcomp


#snpper = robjects.r('function(snpregion,chrm){\n'
#				'library(biomaRt)\n'
#				'snpmart = useMart("snp", dataset="hsapiens_snp")\n'
#				'snprecords = getBM(c("chrom_start","refsnp_id","allele","refsnp_source","validated","mapweight","study_type"),\n'
 #               'filters = c("chromosomal_region","with_validated","chr_name","variation_source"),\n'
  #              'values = list(snpregion,TRUE,chrm,"dbSNP"), mart = snpmart)\n'
	#			'return(snprecords)}')

snpbim = matrix_reader('/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.bim')
doubleup = open('/mnt/lustre/home/cusanovich/snp138Mult.rsIDs.txt')
doublemint = {}
for line in doubleup:
	doublemint[line.strip()] = ''

doubleup.close()

#oldid = snpbim[:,1].copy()
#newid = snpbim[:,1].copy()
attributes = robjects.StrVector(["chrom_start","refsnp_id","allele","refsnp_source","validated","mapweight","study_type"])
nameupdates = open('/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.chr' + chrmreq + '.snpid.updates.txt','w')
#import time
#t0 = time.time()
for item in range(snpbim.shape[0]):
#for item in range(20):
	snp = snpbim[item,]
	chrm = snp[0]
	if str(chrm) != str(chrmreq):
		continue
	oldid = snp[1]
	if item % 1000 == 0:
		print str(item) + ' snps updated...'
		sys.stdout.flush()
	chrm = snp[0]
	snpregion = snp[0] + ':' + snp[3] + ':' + snp[3]
	filters = robjects.Vector([snpregion,True,chrm,"dbSNP"])
	filters.names = ["chromosomal_region","with_validated","chr_name","variation_source"]
	snprecords = biomart.getBM(attributes,filters=filters,mart=snpmart)
	#snprecords = snpper(snp[0] + ':' + snp[3] + ':' + snp[3],snp[0])
	if len(snprecords[0]) == 0:
		print >> nameupdates, oldid + '\tncgi' + snp[1]
		continue
	counter = range(len(snprecords[0]))
	counterupdate = []
	for x,y in enumerate(counter):
		if snprecords[5][y] == 1:
			counterupdate.append(x)
	if len(counterupdate) == 0:
		print >> nameupdates, oldid + '\tmcgi' + snp[1]
		continue
	counter = [counter[x] for x in counterupdate]
	counterupdate = []
	X = snp[4]
	Y = snp[5]
	for x,y in enumerate(counter):
		snpalleles = snprecords[2][y].split('/')
		A = snpalleles[0]
		B = snpalleles[1]
		if A == X:
			if B == Y:
				counterupdate.append(x)
				continue
		if A == Y:
			if B == X:
				counterupdate.append(x)
				continue
		if reverse_complementer(A) == X:
			if reverse_complementer(B) == Y:
				counterupdate.append(x)
				continue
		if reverse_complementer(A) == Y:
			if reverse_complementer(B) == X:
				counterupdate.append(x)
				continue
	if len(counterupdate) == 0:
		print >> nameupdates, oldid + '\tmcgi' + snp[1]
		continue
	counter = [counter[x] for x in counterupdate]
	counterupdate = []
	for x,y in enumerate(counter):
		try:
			h = doublemint[snprecords[1][y]]
		except KeyError:
			counterupdate.append(x)
	if len(counterupdate) == 0:
		print >> nameupdates, oldid + '\tmcgi' + snp[1]
		continue
	counter = [counter[x] for x in counterupdate]
	counterupdate = []
	if len(counter) == 1:
		print >> nameupdates, oldid + '\t' + snprecords[1][counter[0]]
		continue
	for x,y in enumerate(counter):
		if snprecords[6][y] == 'GWAS':
			counterupdate.append(x)
	if len(counterupdate) > 0:
		counter = [counter[x] for x in counterupdate]
		counterupdate = []
	if len(counter) == 1:
		print >> nameupdates, oldid + '\t' + snprecords[1][counter[0]]
		continue
	for x,y in enumerate(counter):
		rsidcount = list(snprecords[1]).count(snprecords[1][y])
		if rsidcount == 1:
			counterupdate.append(x)
		if rsidcount > 1:
			counterupdate.append(min(snprecords[1].index(snprecords[1][y])))
	if len(counterupdate) == 1:
		newid[item] = snprecords[1][counter[counterupdate]]
		continue
	print >> nameupdates, oldid + '\tmcgi' + snp[1]

#print time.time() - t0

#nameupdates = open('/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.snpid.updates.txt','w')
#for update in range(len(oldid)):
#	print >> nameupdates, str(oldid[update]) + '\t' + newid[update]

nameupdates.close()

finisher = open('/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.chr' + chrmreq + '.done','w')
finisher.close()
