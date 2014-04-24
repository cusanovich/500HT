import glob
import numpy
import sys

gwasdir = "/mnt/lustre/home/cusanovich/500HT/MinalSherlock"
gwass = glob.glob(gwasdir + '/*sherlock*')
sensitivity = open("/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.dnase.bed")
sensitivitydic = {}
print "Building sensitivity dictionary..."
for row in sensitivity:
	if 'chrm' in row:
		tissues = row.strip().split()
		continue
	rower = row.strip().split()
	sensitivitydic[rower[3]] = rower

sensitivity.close()
for gwas in gwass:
	currgwas = open(gwas,'r')
	currpheno = gwas.split('/')[-1].split('.')[0]
	outfile = open(gwasdir + '/' + currpheno + '.enrichment.txt.','w')
	print currpheno
	sys.stdout.flush()
	gwasdic = {}
	print "Building GWAS results dictionary..."
	sys.stdout.flush()
	for line in currgwas:
		liner = line.strip().split()
		gwasdic[liner[0]] = float(liner[1])
	currgwas.close()
	tiles = [0.025,0.05,0.1,0.2,0.4,0.8,1.6,3.2,6.4,12.8,25.6,50]
	tiled = [numpy.percentile(gwasdic.values(),x) for x in tiles]
	tiled.append(1.01)
	currsensdic = {}
	currtots = [0]*len(tiled)
	for tissue in tissues[4:]:
		currsensdic[tissue] = [0]*len(tiled)
	print "Calculating enrichments..."
	sys.stdout.flush()
	for snp in gwasdic.keys():
		snpp = gwasdic[snp]
		olaps = sensitivitydic[snp]
		for quant in range(len(tiled)):
			if snpp < tiled[quant]:
				currtots[quant] += 1
			for tissue in currsensdic.keys():
				if olaps[tissues.index(tissue)] != '0':
					currsensdic[tissue][quant] += 1
	print "Writing results..."
	sys.stdout.flush()
	print >> outfile, 'Tissue\t' + '\t'.join([str(x) for x in tiled])
	for tissue in currsensdic.keys():
		baseline = currsensdic[tissue][-1]/currtots[-1]
		print >> outfile, tissue + '\t' + '\t'.join([str((currsensdic[tissue][x]/currtots[x])/baseline) for x in range(len(currsensdic[tissue]))])
	outfile.close()

