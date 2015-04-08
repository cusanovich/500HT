#!/usr/bin/env python

import subprocess
import glob
import time
import sys
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
from myfuncs import ifier

i = sys.argv[1]
linecounter = subprocess.Popen('wc -l /mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.gccor.newcovcor.ordered.chr' + str(i) + '.genes', shell=True, stdout=subprocess.PIPE)
linecount = int(linecounter.communicate()[0].strip().split()[0])
blocks = linecount/100
for j in range(blocks):
	eqtler = 'echo "python /mnt/lustre/home/cusanovich/500HT/Scripts/trans_mapper.py ' + str(i) + ' ' + str(100*j) + ' ' + str(100*(j+1)) + '" | qsub -l h_vmem=2g -V -o ~/dump/ -e ~/dump/ -N "trans.chr' + str(i) + '.block' + str(j) + '"'
	ifier(eqtler)

remained = linecount%100
eqtler = 'echo "python /mnt/lustre/home/cusanovich/500HT/Scripts/trans_mapper.py ' + str(i) + ' ' + str(100*blocks) + ' ' + str(100*(blocks) + remained) + '" | qsub -l h_vmem=2g -V -o ~/dump/ -e ~/dump/ -N "trans.chr' + str(i) + '.block' + str(blocks) + '"'
ifier(eqtler)

while len(glob.glob('/mnt/lustre/home/cusanovich/500HT/ByChr/chr' + str(i) + '.*.done')) < (blocks + 1):
	time.sleep(300)

cleanup = "rm /mnt/lustre/home/cusanovich/500HT/ByChr/chr" + str(i) + ".*.done"
ifier(cleanup)

masterer = ('cat /mnt/lustre/home/cusanovich/500HT/ByChr/chr' + str(i) +
	'.*.newcovcor.trans.pvals.txt > /mnt/lustre/home/cusanovich/500HT/ByChr/chr' + str(i) +
	'.newcovcor.trans.pvals.txt; cat /mnt/lustre/home/cusanovich/500HT/ByChr/chr' + str(i) +
	'.*.trans.sig.txt >  /mnt/lustre/home/cusanovich/500HT/ByChr/chr' + str(i) +
	'.trans.sig.txt; cat /mnt/lustre/home/cusanovich/500HT/ByChr/chr' + str(i) +
	'.*.newcovcor.sherlock.txt >  /mnt/lustre/home/cusanovich/500HT/ByChr/chr' + str(i) +
	'.newcovcor.sherlock.txt')
ifier(masterer)

cleanup = 'rm /mnt/lustre/home/cusanovich/500HT/Imputed1415/chr' + i + '.fam; rm /mnt/lustre/home/cusanovich/500HT/Imputed1415/chr' + i + '.bed; rm /mnt/lustre/home/cusanovich/500HT/Imputed1415/chr' + i + '.bim; rm /mnt/lustre/home/cusanovich/500HT/ByChr/chr' + i + '.block*'
ifier(cleanup)

zipper = 'gzip /mnt/lustre/home/cusanovich/500HT/ByChr/chr' + str(i) + '.newcovcor.trans.pvals.txt'
ifier(zipper)
