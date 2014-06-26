import sys
import os
import numpy
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
from DarrenTools import ifier, matrix_reader
genodir = '/mnt/lustre/home/cusanovich/500HT/Imputed1415/'
os.chdir(genodir)

plinker = 'plink --bfile 500HT/Imputed1415/hutt.imputed.rename --recodeA'
ifier(plinker)

genomatix = matrix_reader('/mnt/lustre/home/cusanovich/plink.raw',sep=" ")
properg = genomatix.T[6:,1:]
snps = [x.split('_')[0] for x in genomatix.T[6:,0]]
findivs = genomatix.T[1,1:]

numpy.savetxt('/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.rename.genos.bimbam', properg, fmt = '%s', delimiter = '\t', newline = '\n')

snper = open('/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.rename.snps.bimbam','w')
print >> snper, '\n'.join(snps)
snper.close()

finer = open('/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.rename.findivs.bimbam','w')
print >> finer, '\n'.join(findivs)
finer.close()

