#!/usr/bin/env python
import subprocess
import time
import glob


#step 1 - chrom beds
#step 2 - chrom raws
#step 3 - chrom bim (snpid + position)
#step 4 - load raw
#step 5 - transpose raw
#step 6 - add chrom, start, end, snpid, "0", "."
#step 7 - save table
#step 8 - compress with bgzip
#step 9 - index with tabix

def ifier(commander):
        ify = subprocess.Popen(commander,shell=True)
        ify.wait()

genodir = '/mnt/lustre/home/cusanovich/500HT/Imputed1415/'

# print 'Creating raw files...'
# for j in range(1,23):
# 	plinker = 'echo "plink --noweb --nonfounders --maf 0.05 --geno 0.05 --bfile ' + genodir + 'imputed_cgi --chr ' + str(j) + ' --make-bed --out ' + genodir + 'ByChr/hutt.imputed.chr' + str(j) + '; plink --bfile ' + genodir + 'ByChr/hutt.imputed.chr' + str(j) + ' --recodeA --out ' + genodir + 'ByChr/hutt.imputed.chr' + str(j) + '; touch ' + genodir + 'ByChr/chr' + str(j) + '.done" | qsub -l h_vmem=2g -o ~/dump/ -e ~/dump/'
# 	ifier(plinker)

# while len(glob.glob(genodir + 'ByChr/*.done')) < 22:
# 	time.sleep(5)

# cleanup = "rm " + genodir + "ByChr/*.done"
# ifier(cleanup)

# print 'Creating bed files...'
# for j in range(1,23):
# 	converter = 'echo "python /mnt/lustre/home/cusanovich/500HT/Scripts/raw2txt.py ' + str(j) + ' ' + genodir + '" | qsub -l h_vmem=8g -o ~/dump/ -e ~/dump/'
# 	ifier(converter)

# while len(glob.glob(genodir + 'ByChr/*.done')) < 22:
# 	time.sleep(5)

# ifier(cleanup)

print 'Creating tabix files...'
for j in range(1,23):
	tabixer = 'echo "bgzip ' + genodir + 'ByChr/hutt.imputed.chr' + str(j) + '.txt; tabix -p bed ' + genodir + 'ByChr/hutt.imputed.chr' + str(j) + '.txt.gz; touch ' + genodir + 'ByChr/' + str(j) + '.done" | qsub -l h_vmem=2g -o ~/dump/ -e ~/dump/ -wd `pwd`'
	ifier(tabixer)

while len(glob.glob(genodir + 'ByChr/*.done')) < 22:
	time.sleep(5)

print 'Cleaning up a bit...'
ifier(cleanup)
for j in range(1,23):
	for k in ['.bed','.bim','.fam','.log','.nof','.raw']:
		cleaner = 'rm ' + genodir + 'ByChr/hutt.imputed.chr' + str(j) + k
		ifier(cleaner)

print 'Work complete.'