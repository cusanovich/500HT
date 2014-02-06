#!/usr/bin/env python

import sys
import subprocess
import glob
import time

hmdir = '/mnt/lustre/home/cusanovich/'
pc = str(sys.argv[1])

def ifier(commander):
	ify = subprocess.Popen(commander,shell=True)
	ify.wait()

#for pc in range(21):
#for pc in [4] + range(10,21):
print pc
print "Regressing PCs..."
for block in range(15):
	gemmer = 'echo "python /mnt/lustre/home/cusanovich/500HT/Scripts/gemma.py ' + str(block) + ' ' + str(pc) + '" | qsub -l h_vmem=2g -o ~/dump/ -e ~/dump/'
	ifier(gemmer)
while len(glob.glob(hmdir + '500HT/output/*.PC' + str(pc) + '.residE.txt')) < 14055:
	time.sleep(60)
print "Combining regressions..."
for i in range(15):
	j = str(i)
	if i < 10:
		j = '0' + str(i)
	paster = 'paste ' + hmdir + '500HT/output/*.Block' + j + '.PC' + str(pc) + '.residE.txt > ' + hmdir + '500HT/output/' + j + '.PC' + str(pc) + '.qqnorm.500ht.bimbam.fixed'
	ifier(paster)
print "Cleaning up..."
paster = 'paste ' + hmdir + '500HT/output/*.PC' + str(pc) + '.qqnorm.500ht.bimbam.fixed > ' + hmdir + '500HT/Exprs/qqnorm.500ht.bimbam.PC' + str(pc) + '.fixed; find ' + hmdir + '500HT/output/ -name "*.PC' + str(pc) + '.*" -exec rm {} \;'
ifier(paster)