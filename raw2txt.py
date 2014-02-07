#!/usr/bin/env python
import sys
import subprocess
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
import numpy

def ifier(commander):
	ify = subprocess.Popen(commander,shell=True)
	ify.wait()

chrm = str(sys.argv[1])
genodir = sys.argv[2]

genor = genodir + 'ByChr/hutt.3chip.chr' + chrm + '.genos.txt'
anoter = genodir + 'ByChr/hutt.3chip.chr' + chrm + '.annot.txt'

print "Loading annotations..."
sys.stdout.flush()
bim = open(genodir + 'ByChr/hutt.3chip.chr' + chrm + '.bim','r')
bimer = bim.readlines()
snpids = [x.strip().split()[1] for x in bimer]
snppos = [x.strip().split()[3] for x in bimer]

print "Loading genotypes..."
sys.stdout.flush()
linecounter = subprocess.Popen('wc -l ' + genodir + 'ByChr/hutt.3chip.chr' + chrm + '.raw', shell=True, stdout=subprocess.PIPE)
linecount = int(linecounter.communicate()[0].strip().split()[0])
columncounter = subprocess.Popen('awk -F" " \'{print NF;exit}\' ' + genodir + 'ByChr/hutt.3chip.chr' + chrm + '.raw', shell=True, stdout=subprocess.PIPE)
columncount = int(columncounter.communicate()[0].strip().split()[0])
raws = numpy.zeros((linecount,columncount),dtype='|S2')
rawin = open(genodir + 'ByChr/hutt.3chip.chr' + chrm + '.raw','r')
for i,line in enumerate(rawin):
	raws[i,:] = line.strip().split()

#raws = numpy.loadtxt(genodir + 'ByChr/hutt.imputed.chr' + chrm + '.raw',dtype='str')

print "Transposing genotypes..."
sys.stdout.flush()
genorfile = open(genor,'w')
for line in range(6,raws.shape[1]):
	print >> genorfile, "\t".join(list(raws[:,line]))

genorfile.close()

#traws = numpy.transpose(raws)
#trawsu = traws[6:traws.shape[0],1:traws.shape[1]]
#numpy.savetxt(genor,trawsu,delimiter="\t",fmt='%s')

print "Annotating SNPs..."
sys.stdout.flush()
traws1 = numpy.column_stack((['.']*len(snpids),[0]*len(snpids)))
traws2 = numpy.column_stack((snpids,traws1))
traws3 = numpy.column_stack((snppos,traws2))
traws4 = numpy.column_stack(([str(int(x)-1) for x in snppos],traws3))
traws5 = numpy.column_stack((['chr' + str(chrm)]*len(snpids),traws4))
numpy.savetxt(anoter,traws5,delimiter="\t",fmt='%s')

print "Finalizing files..."
sys.stdout.flush()
paster = '/bin/bash -c "paste <(cat ' + anoter + ') <(cat ' + genor + ') > ' + genodir + 'ByChr/hutt.3chip.chr' + chrm + '.txt; rm ' + anoter + '; rm ' + genor + '"'
ifier(paster)

doner = open(genodir + 'ByChr/chr' + chrm + '.done','w')
doner.close()