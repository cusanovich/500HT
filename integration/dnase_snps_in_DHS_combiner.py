import sys
sys.path.insert(0,'/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/pybedtools-0.6.2-py2.6-linux-x86_64.egg/pybedtools')
import pybedtools
from pybedtools import BedTool
import glob
import subprocess

files = glob.glob("/mnt/lustre/home/cusanovich/500HT/GenomeAnnotations/*.bed")
firstbed = BedTool(files[0])
snpbed = BedTool("/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.snp_olaping_dhs_coords.bed")
snpbed.intersect(firstbed,c=True).saveas("/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.dnase2.bed")
titler = 'echo "chrm\tstart\tend\trsID\t' + files[0].split('/')[-1].split('.')[0] + '" | cat - /mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.dnase2.bed > /mnt/lustre/home/cusanovich/500HT/Imputed1415/last_tissue.bed'
titling = subprocess.Popen(titler, shell=True)
titling.wait()
for x,fname in enumerate(files[1:]):
	print str(x)
	sys.stdout.flush()
	tissue = fname.split("/")[-1].split(".")[0]
	print tissue
	currbed = BedTool(fname)
	snpbed = BedTool("/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.snp_olaping_dhs_coords.bed")
	snpbed.intersect(currbed,c=True).saveas("/mnt/lustre/home/cusanovich/500HT/Imputed1415/curr.bed")
	newcol = 'echo "' + tissue + '" > /mnt/lustre/home/cusanovich/500HT/Imputed1415/curr.count.txt; cut -f5 /mnt/lustre/home/cusanovich/500HT/Imputed1415/curr.bed >> /mnt/lustre/home/cusanovich/500HT/Imputed1415/curr.count.txt; paste /mnt/lustre/home/cusanovich/500HT/Imputed1415/last_tissue.bed /mnt/lustre/home/cusanovich/500HT/Imputed1415/curr.count.txt > /mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.dnase2.bed; cp /mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.dnase2.bed /mnt/lustre/home/cusanovich/500HT/Imputed1415/last_tissue.bed'
	print newcol
	newcoling = subprocess.Popen(newcol, shell=True)
	newcoling.wait()
	cleaner = 'rm /mnt/lustre/home/cusanovich/500HT/Imputed1415/curr.count.txt; rm /mnt/lustre/home/cusanovich/500HT/Imputed1415/curr.bed'
	cleaning = subprocess.Popen(cleaner, shell=True)
	cleaning.wait()

cleaner = 'rm /mnt/lustre/home/cusanovich/500HT/Imputed1415/last_tissue.bed'
cleaning = subprocess.Popen(cleaner, shell=True)
cleaning.wait()
