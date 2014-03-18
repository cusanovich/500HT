#!/usr/bin/env python

#import subprocess
#import glob
#import time
import sys
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
from myfuncs import ifier

for i in range(22,0,-1):
	eqtler = 'nohup python /mnt/lustre/home/cusanovich/500HT/Scripts/trans_eqtl_driver.py ' + str(i) + ' &'
	ifier(eqtler)
	#eqtler = 'echo ' + str(i)
	#ifier(eqtler)
