#!/usr/bin/env python

import subprocess

def ifier(commander):
        ify = subprocess.Popen(commander,shell=True)
        ify.wait()

for i in range(1,9):
        eqtler = 'echo "python /mnt/lustre/home/cusanovich/500HT/Scripts/rs_grabber.py ' + str(i) + '" | qsub -l h_vmem=2g -V -o ~/dump/ -e ~/dump/ -N "rsIDs.chr' + str(i) + '"'
        ifier(eqtler)
