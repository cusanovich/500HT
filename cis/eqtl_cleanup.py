#!/usr/bin/env python

import subprocess
import glob
import time
import sys

def ifier(commander):
	ify = subprocess.Popen(commander,shell=True)
	ify.wait()

for pcs in range(0,101):
	print pcs
	if len(glob.glob('/mnt/lustre/home/cusanovich/500HT/ByChr/*.PC' + str(pcs) + '.bonferroni.done')) == 22:
		cleanup = "rm /mnt/lustre/home/cusanovich/500HT/ByChr/*.PC" + str(pcs) + ".bonferroni.done"
		ifier(cleanup)
		masterer = 'cat /mnt/lustre/home/cusanovich/500HT/ByChr/*.PC' + str(pcs) + '.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.eqtls.txt | sort -k1,1 -k2,2 > /mnt/lustre/home/cusanovich/500HT/eQTLs/master.PC' + str(pcs) + '.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.eqtls.txt; cat /mnt/lustre/home/cusanovich/500HT/ByChr/*.PC' + str(pcs) + '.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.chosen.txt | sort -k1,1 >  /mnt/lustre/home/cusanovich/500HT/eQTLs/master.PC' + str(pcs) + '.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.chosen.txt'
		ifier(masterer)
