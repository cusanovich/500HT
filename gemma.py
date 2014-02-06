#!/usr/bin/env python

import sys
import subprocess

hmdir = '/mnt/lustre/home/cusanovich/'
multiplier = sys.argv[1]
pcs = sys.argv[2]
block = str(multiplier)
if int(multiplier) < 10:
	block = '0' + str(multiplier)
filetag = '.Block' + block + '.PC' + str(pcs)

def ifier(commander):
	ify = subprocess.Popen(commander,shell=True)
	ify.wait()

def gemming(i):
	lister = str(i)
	if len(lister) < 5:
		lister = '0'*(5-len(lister)) + lister
	gemmer = ('echo ' + str(i) + '; cd ' + hmdir + '500HT/; ' + hmdir +
			  'Programs/gemma0.94 -p ' + hmdir +
			  '500HT/qqnorm.500ht.bimbam -k ' + hmdir +
			  '500HT/addSNP.coef.3671.square -c ' + hmdir +
			  '500HT/Exprs/qqnorm.500ht.pc' + str(pcs) + ' -lmm 5 -n ' +
			  str(i) + ' -o ' + lister + filetag + '; rm ' + hmdir +
			  '500HT/output/' + lister + filetag + '.log.txt' + '; rm ' +
			  hmdir + '500HT/output/' + lister + filetag + '.residU.txt')
	ifier(gemmer)

for i in range((int(multiplier)*1000+1),(int(multiplier)*1000+1001)):
	if i > 14055:
		continue
	gemming(i)