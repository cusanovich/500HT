import glob
import os

outfile = open('./data/RNAseq_scripts/readsbylane.matrix','w')
print >> outfile, 'Total\tDemultied\tSaved\tMapped\tJuncMapped\tSaveMapped\tJuncsaved'
for x in range(1,12):
	currfc = 'FlowCell' + str(x)
	currfcid = 'fc' + str(x)
	if x < 10:
		currfcid = 'fc0' + str(x)
	findivs = os.listdir('./data/500HTRNAseq/' + currfc)
	for findiv in findivs:
		counts = glob.glob('./data/500HTRNAseq/' + currfc + '/' + findiv + '/*.counts.txt')
		for count in counts:
			currcount = open(count,'r').readlines()
			ider = findiv + '.' + count.split('/')[-1].split('.')[0] + '.' + currfcid
			tot = str(int(currcount[0].strip().split()[1]) + int(currcount[5].strip().split()[1]))
			totd = currcount[0].strip().split()[1]
			tots = currcount[5].strip().split()[1]
			mapped = currcount[1].strip().split()[1]
			juncmapped = currcount[3].strip().split()[1]
			savemapped = currcount[6].strip().split()[1]
			juncsaved = currcount[7].strip().split()[1]
			currline = [ider, tot, totd, tots, mapped, juncmapped, savemapped, juncsaved]
			print >> outfile, '\t'.join(currline)

outfile.close()
