#!/usr/bin/env python
import os
import sys
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
from myfuncs import genereader
import subprocess
import numpy
import pysam
import gzip
import rpy2.robjects as robjects

# def ifier(commander):
# 	ify = subprocess.Popen(commander,shell=True)
# 	ify.wait()

# def matrix_reader(matrix_file,sep="\t",dtype='|S20'):
# 	linecounter = subprocess.Popen('wc -l ' + matrix_file, shell=True, stdout=subprocess.PIPE)
# 	linecount = int(linecounter.communicate()[0].strip().split()[0])
# 	columncounter = subprocess.Popen('awk -F"' + sep + '" \'{print NF;exit}\' ' + matrix_file, shell=True, stdout=subprocess.PIPE)
# 	columncount = int(columncounter.communicate()[0].strip().split()[0])
# 	raws = numpy.zeros((linecount,columncount),dtype=dtype)
# 	rawin = open(matrix_file,'r')
# 	for i,line in enumerate(rawin):
# 		raws[i,:] = line.strip().split()
# 	rawin.close()
# 	return raws

corr = robjects.r('function(master,target){\n'
				 'x = as.numeric(unlist(exprs))\n'
				 'y = as.numeric(unlist(genos))\n'
				 'corr = cor(x,y,method="spearman")\n'
				 'return(corr)}')

pcs = 60
something to read in list of tfs here...
grepdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/Union/10kb/'
genelists = os.listdir(grepdir)
genelist = sorted(set([x.split('.')[0] for x in genelists]))
targets = {}
mastergene = {}
for tf in genelist:
	try:
		targets[tf] = genereader(grepdir+gene+'.DEandBound.txt')
		mastergene[tf] = genereader(grepdir+gene+'.BindingTF.txt')
	except IOError:
		pass

expressed = matrix_reader('/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.covcor.ordered.bimbam',dtype='float')
genes = matrix_reader('/mnt/lustre/home/cusanovich/500HT/genenames.500ht.txt')
findivs = matrix_reader('/mnt/lustre/home/cusanovich/500HT/findivs.500ht.ordered.txt')
masterind = {}
targetinds = {}
for gene in mastergene.keys():
	masterind[gene] = genes.index(mastergene[gene])
	targetinds[gene] = [genes.index(x) for x in targets[gene]]

totalcor = []
for pc in range(pcs):
	pcmat = matrix_reader('/mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.gccor.covcor.ordered.pc' + str(pc),dtype='float')
	mod1 = pcmat
	mod2 = mod1.T
	Y = expressed.T
	if pc == 0:
		Yfit = Y
	else:
		W = pcmat[:,1:pc]
		mods = mod2.dot(mod1)
		inversemods = numpy.linalg.inv(mods)
		gammahat = Y.dot(mod1).dot(inversemods)[:,1:mods.shape[1]]
		Yfit = Y - gammahat.dot(W.T)
	tfcor = []
	for tf in targetinds.keys():
		master = robjects.FloatVector(Yfit[masterind[tf],:])
		targetcor = []
		for row in targetinds[tf]:
			targetcor.append(corr(master,robjects.FloatVector(Yfit[row,:])))
		tfcor.append(numpy.median(targetcor))
	totalcor.append(numpy.median(targetcor))

print max(totalcor)
print totalcor.index(max(totalcor))