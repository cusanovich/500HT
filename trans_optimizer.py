#!/usr/bin/env python
import os
import sys
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
from myfuncs import genereader,matrix_reader
import subprocess
import numpy
import pysam
import gzip
import rpy2.robjects as robjects

corr = robjects.r('function(master,target){\n'
				 'corr = cor.test(master,target,method="spearman",exact=FALSE)\n'
				 'return(c(abs(corr$estimate),corr$p.value))}')

padjcounter = robjects.r('function(ps){\n'
				 'padjs = p.adjust(ps,method="BH")\n'
				 'counted = length(which(padjs < 0.05))\n'
				 'return(counted)}')

def bonfcounter(ps):
	adjuster = len(ps)
	counted = 0
	for value in ps:
		if value*adjuster < 0.05:
			counted += 1
	return(counted)

pcs = 60
mastergene = {}
tflist = open('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt','r')
for line in tflist:
	if '#' in line:
		continue
	try:
		test = mastergene[line.strip().split()[0]]
	except KeyError:
		mastergene[line.strip().split()[0]] = line.strip().split()[2]

grepdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Bindings/'
genelists = os.listdir(grepdir)
genelist = sorted(set([x.split('.')[0] for x in genelists]))
targets = {}
for tf in genelist:
	try:
		targets[tf] = genereader(grepdir+tf+'.DEandBound.txt')
	except IOError:
		pass

expressed = matrix_reader('/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.covcor.ordered.bimbam',dtype='float')
genes = open('/mnt/lustre/home/cusanovich/500HT/genenames.500ht.txt','r').readlines()
genes = [x.strip() for x in genes]
genedic = {}
for g,k in enumerate(genes):
	genedic[k] = g

masterind = {}
targetinds = {}
for gene in targets.keys():
	try:
		masterind[gene] = genedic[mastergene[gene]]
	except KeyError:
		continue
	targetinds[gene] = []
	for subgene in targets[gene]:
		try:
			targetinds[gene].append(genedic[subgene])
		except KeyError:
			continue

totalcor = []
totalps = []
for pc in range(pcs):
	print str(pc) + ' PCs are being removed from the expression data...'
	sys.stdout.flush()
	pcmat = matrix_reader('/mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.gccor.covcor.ordered.pc' + str(pc),dtype='float')
	Y = expressed.T
	if pc == 0:
		Yfit = Y
	else:
		mod1 = pcmat
		mod2 = mod1.T
		W = pcmat[:,1:(pc+1)]
		mods = mod2.dot(mod1)
		inversemods = numpy.linalg.inv(mods)
		gammahat = Y.dot(mod1).dot(inversemods)[:,1:mods.shape[1]]
		Yfit = Y - gammahat.dot(W.T)
	tfcor = []
	tfps = []
	for tf in targetinds.keys():
		master = robjects.FloatVector(Yfit[masterind[tf],:])
		targetcor = []
		for row in targetinds[tf]:
			corring = corr(master,robjects.FloatVector(Yfit[row,:]))
			targetcor.append(corring[0])
			tfps.append(corring[1])
		tfcor.append(numpy.median(targetcor))
	totalcor.append(numpy.median(targetcor))
	print totalcor[-1]
	totalps.append(bonfcounter(tfps))
	print totalps[-1]
	print min(tfps)

outfile = open('/mnt/lustre/home/cusanovich/500HT/trans_maximizer.txt','w')
for z in range(len(totalcor)):
	print >> outfile, str(totalcor[z]) + '\t' + str(totalps[z])

outfile.close()

print max(totalcor)
print totalcor.index(max(totalcor))
print max(totalps)
print totalps.index(max(totalps))
