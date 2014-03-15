#!/usr/bin/env python

import subprocess
import glob
import time
import sys
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
from myfuncs import ifier

for i in range(22,0,-1):
	eqtler = 'echo "python /mnt/lustre/home/cusanovich/500HT/Scripts/trans_mapper.py ' + str(i) + '" | qsub -l h_vmem=2g, -V -o ~/dump/ -e ~/dump/ -N "trans_eQTLs.chr' + str(i) + '"'
	ifier(eqtler)

#while len(glob.glob('/mnt/lustre/home/cusanovich/500HT/ByChr/*.PC' + str(pcs) + '.bonferroni.done')) < 22:
#	time.sleep(300)

#cleanup = "rm /mnt/lustre/home/cusanovich/500HT/ByChr/*.PC" + str(pcs) + ".bonferroni.done"
#ifier(cleanup)

#masterer = 'cat /mnt/lustre/home/cusanovich/500HT/ByChr/*.PC' + str(pcs) + '.imputed.150kb.bonferroni.gccor.covcor.regressPCs.gemma.eqtls.txt | sort -k1,1 -k2,2 > /mnt/lustre/home/cusanovich/500HT/eQTLs/master.PC' + str(pcs) + '.imputed.150kb.bonferroni.gccor.covcor.regressPCs.gemma.eqtls.txt; cat /mnt/lustre/home/cusanovich/500HT/ByChr/*.PC' + str(pcs) + '.imputed.150kb.bonferroni.gccor.covcor.regressPCs.gemma.chosen.txt | sort -k1,1 >  /mnt/lustre/home/cusanovich/500HT/eQTLs/master.PC' + str(pcs) + '.imputed.150kb.bonferroni.gccor.covcor.regressPCs.gemma.chosen.txt'
#ifier(masterer)
