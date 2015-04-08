infile = open("/mnt/lustre/home/cusanovich/500HT/eQTLs/remaster.imputed.1Mb.bonferroni.gccor.newcovcor.sherlock.txt","r")
outfile = open("/mnt/lustre/home/cusanovich/500HT/eQTLs/remaster.imputed.1Mb.bonferroni.gccor.newcovcor.sherlock.filter.txt","w")
for line in infile:
	liner = line.strip().split()
	if liner[3] == '0' and float(liner[2]) < 1.0e-05:
		print >> outfile, "\t".join(liner)
	if liner[3] == '1' and float(liner[2]) < 1.0e-03:
		print >> outfile, "\t".join(liner)

infile.close()
outfile.close()
