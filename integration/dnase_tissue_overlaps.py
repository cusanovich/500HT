import sys
sys.path.insert(0,'/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/pybedtools-0.6.2-py2.6-linux-x86_64.egg/pybedtools')
import pybedtools
from pybedtools import BedTool
import glob
import subprocess

files = glob.glob("/mnt/lustre/home/cusanovich/500HT/GenomeAnnotations/*.bed")
firstbed = BedTool(files[0])
masterbed = BedTool("/mnt/lustre/home/cusanovich/500HT/master_dhs_merged.bed")
masterbed.intersect(firstbed,c=True).saveas("/mnt/lustre/home/cusanovich/500HT/master.dnase.overlap.bed")
titler = 'echo "chrm\tstart\tend\t' + files[0].split('/')[-1].split('.')[0] + '" | cat - /mnt/lustre/home/cusanovich/500HT/master.dnase.overlap.bed > /mnt/lustre/home/cusanovich/500HT/last_tissue.bed'
titling = subprocess.Popen(titler, shell=True)
titling.wait()
for x,fname in enumerate(files[1:]):
	print str(x)
	sys.stdout.flush()
	tissue = fname.split("/")[-1].split(".")[0]
	print tissue
	currbed = BedTool(fname)
	masterbed = BedTool("/mnt/lustre/home/cusanovich/500HT/master_dhs_merged.bed")
	masterbed.intersect(currbed,c=True).saveas("/mnt/lustre/home/cusanovich/500HT/curr.bed")
	newcol = 'echo "' + tissue + '" > /mnt/lustre/home/cusanovich/500HT/curr.count.txt; cut -f4 /mnt/lustre/home/cusanovich/500HT/curr.bed >> /mnt/lustre/home/cusanovich/500HT/curr.count.txt; paste /mnt/lustre/home/cusanovich/500HT/last_tissue.bed /mnt/lustre/home/cusanovich/500HT/curr.count.txt > /mnt/lustre/home/cusanovich/500HT/master.dnase.overlap.bed; cp /mnt/lustre/home/cusanovich/500HT/master.dnase.overlap.bed /mnt/lustre/home/cusanovich/500HT/last_tissue.bed'
	print newcol
	newcoling = subprocess.Popen(newcol, shell=True)
	newcoling.wait()
	cleaner = 'rm /mnt/lustre/home/cusanovich/500HT/curr.count.txt; rm /mnt/lustre/home/cusanovich/500HT/curr.bed'
	cleaning = subprocess.Popen(cleaner, shell=True)
	cleaning.wait()

cleaner = 'rm /mnt/lustre/home/cusanovich/500HT/last_tissue.bed'
cleaning = subprocess.Popen(cleaner, shell=True)
cleaning.wait()
