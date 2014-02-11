#!/usr/bin/env python

print "Loading SNP names..."
#genos = open('/mnt/lustre/home/cusanovich/500HT/Imputed1415/imputed_cgi.bim','r')
genos = open('/mnt/lustre/home/cusanovich/oldhome_nb/Hutterite_Heritability/gemma/Genotypes/hutt.3chip.hg19remap.bim','r')

genonamed = genos.readlines()
genonamer = [x.strip().split()[0] for x in genonamed]
genonamed = [x.strip().split()[1] for x in genonamed]
genonames = {}
genonamers = {}
for x, j in enumerate(genonamed):
	genonames[j] = x
	genonamers[j] = 'chr' + str(genonamer[x])

genos.close()

print "Loading gene names..."
naming = open('/mnt/lustre/home/cusanovich/500HT/genenames.500ht.txt','r')
exprnamed = naming.readlines()
exprnamed = [x.strip().split()[0] for x in exprnamed]
exprnames = {}
for x, j in enumerate(exprnamed):
	exprnames[j] = x

naming.close()

chrmexprnames = {}
for chrm in range(1,23):
	currchrm = 'chr' + str(chrm)
	naming = open('/mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.3chip_order.' + currchrm + '.genes','r')
	exprnamed = naming.readlines()
	exprnamed = [x.strip().split()[0] for x in exprnamed]
	chrmexprnames[currchrm] = {}
	for x, j in enumerate(exprnamed):
		chrmexprnames[currchrm][j] = x
	naming.close()

outfile = open('/mnt/lustre/home/cusanovich/500HT/hutt.3chip.mastercols.txt','w')
chroutfile = open('/mnt/lustre/home/cusanovich/500HT/hutt.3chip.chrmspecific.mastercols.txt','w')

print "Recording columns..."
overlaps = open('/mnt/lustre/home/cusanovich/500HT/hutt.3chip.overlap.txt','r')
for line in overlaps:
	liner = line.strip().split()
	chrm = genonamers[liner[1]]
	print >> outfile, liner[0] + "\t" + liner[1] + "\t" + str(exprnames[liner[0]]) + "\t" + str(genonames[liner[1]]) + "\t" + genonamers[liner[1]]
	print >> chroutfile, liner[0] + "\t" + liner[1] + "\t" + str(chrmexprnames[chrm][liner[0]]) + "\t" + str(genonames[liner[1]]) + "\t" + genonamers[liner[1]]

overlaps.close()
outfile.close()
chroutfile.close()
