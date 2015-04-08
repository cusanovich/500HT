import random
from random import shuffle
import sys

seeder = sys.argv[1]


realgenos = open("/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.rename.genos.bimbam","r")
realsnps = open("/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.rename.snps.bimbam","r")
currperm = open("/mnt/lustre/home/cusanovich/500HT/Imputed1415/perms/perm." + seeder + ".genos.bimbam","w")
randind = range(1415)
random.seed(seeder)
random.shuffle(randind)

for line in realgenos:
	liner = line.strip().split()
	newgenos = [liner[ind] for ind in randind]
	currsnp = realsnps.readline().strip()
	mins = liner.count("2")
	majs = liner.count("0")
	hets = liner.count("1")
	maf = round(float(mins*2 + hets)/(2*(mins + majs + hets)),3)
	newline = [currsnp,"A","G",str(maf)] + newgenos
	print >> currperm, ",".join(newline)

realgenos.close()
realsnps.close()
currperm.close()
