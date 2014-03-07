#!/usr/bin/env python

import os
import glob
import subprocess

def ifier(commander):
	ify = subprocess.Popen(commander,shell=True)
	ify.wait()

for pc in range(0,41):
	if len(glob.glob("/mnt/lustre/home/cusanovich/500HT/ByChr/*.PC" + str(pc) + ".bonferroni.done")) == 0:
		continue
	for chrm in range(22,0,-1):
		if os.path.isfile("/mnt/lustre/home/cusanovich/500HT/ByChr/chr" + str(chrm) + ".PC" + str(pc) + ".bonferroni.done"):
			continue
		eqtler = 'echo "python /mnt/lustre/home/cusanovich/500HT/Scripts/alt_gemma_noplink_eqtl_mapper.py chr' + str(chrm) + ' ' + str(pc) + '" | qsub -l h_vmem=2g,bigio=4 -V -o ~/dump/ -e ~/dump/ -N "eQTLs.chr' + str(chrm) + '.PC' + str(pc) + '"'
		ifier(eqtler)
