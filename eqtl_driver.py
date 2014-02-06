#!/usr/bin/env python

import subprocess
import glob
import time
import sys

if len(sys.argv) != 2:
	print "Usage = python eqtl_driver.py pc"
	sys.exit()

pcs = sys.argv[1]

def ifier(commander):
	ify = subprocess.Popen(commander,shell=True)
	ify.wait()

for i in range(22,0,-1):
	eqtler = 'echo "python /mnt/lustre/home/cusanovich/500HT/Scripts/alt_gemma_eqtl_mapper.py chr' + str(i) + ' ' + str(pcs) + '" | qsub -l h_vmem=6g -o ~/dump/ -e ~/dump/ -N "eQTLs.chr' + str(i) + '.PC' + str(pcs) + '"'
	ifier(eqtler)

while len(glob.glob('/mnt/lustre/home/cusanovich/500HT/ByChr/*.PC' + str(pcs) + '.done')) < 22:
	time.sleep(300)

cleanup = "rm /mnt/lustre/home/cusanovich/500HT/ByChr/*.PC" + str(pcs) + ".done"
ifier(cleanup)

masterer = 'cat /mnt/lustre/home/cusanovich/500HT/ByChr/*.PC' + str(pcs) + '.3chip.gemma.eqtls.txt | sort -k1,1 -k2,2 > /mnt/lustre/home/cusanovich/500HT/eQTLs/master.PC' + str(pcs) + '.3chip.gemma.eqtls.txt; cat /mnt/lustre/home/cusanovich/500HT/ByChr/*.PC' + str(pcs) + '.3chip.gemma.chosen.txt | sort -k1,1 >  /mnt/lustre/home/cusanovich/500HT/eQTLs/master.PC' + str(pcs) + '.3chip.gemma.chosen.txt'
ifier(masterer)
