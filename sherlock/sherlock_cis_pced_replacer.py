#!/usr/bin/env python
import os
import sys
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
from DarrenTools import matrix_reader,ifier
import subprocess

#Script to replace uncorrected cis p-values with p-values from PC regressed analysis

badsnps = ["mcgi983490","rs11900082","mcgi2090342","rs7627112","mcgi2118071","rs10865809","mcgi2880951","rs6824720","mcgi3063830","rs73177688","rs2291782","mcgi3913760","rs1990908","mcgi4170170","ncgi4677272","ncgi4677275","mcgi6158705","rs6982255","mcgi6176005","rs2410660","mcgi6384783","rs1431659","mcgi7704523","rs10884261","rs1622217","mcgi7748792","mcgi8038165","rs1532905","mcgi8818315","rs6539584","mcgi9340877","rs9574174","mcgi9598934","rs1168551","mcgi12307275","rs131996"]

cisdic = {}
cisfile = open('/mnt/lustre/home/cusanovich/500HT/eQTLs/master.PC62.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.eqtls.txt','r')
for line in cisfile:
	liner = line.strip().split()
	cisdic[(liner[0],liner[1])] = liner[2]

cisfile.close()
eqtlfile = open('/mnt/lustre/home/cusanovich/500HT/eQTLs/remaster.imputed.1Mb.bonferroni.gccor.newcovcor.sherlock.txt','r')
newfile = open('/mnt/lustre/home/cusanovich/500HT/eQTLs/remaster.imputed.1Mb.bonferroni.gccor.newcovcor.sherlock.cisPCPs.filter.txt','w')

for line in eqtlfile:
	liner = line.strip().split()
	if liner[3] == '0' and float(liner[2]) <= 1e-05:
		print >> newfile, "\t".join(liner)
		continue
	if liner[3] == '1' and liner[1] not in badsnps:
		try:
			#replacer = subprocess.Popen('grep ' + liner[0] + ' /mnt/lustre/home/cusanovich/500HT/eQTLs/master.PC62.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.eqtls.txt | grep ' + liner[1],stdout=subprocess.PIPE,shell=True)
			#newp = replacer.communicate()[0].strip().split()[2]
			#liner[2] = newp
			liner[2] = cisdic[(liner[0],liner[1])]
			if float(liner[2]) <= 1e-03:
				print >> newfile, "\t".join(liner)
		except IndexError:
			print liner
			break

eqtlfile.close()
newfile.close()
