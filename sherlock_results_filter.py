import os
#indir = '/mnt/lustre/home/cusanovich/500HT/MinalSherlock/'
indir = '/mnt/lustre/home/cusanovich/ForMinal/'
infiles = os.listdir(indir)
for infilename in infiles:
	if 'results.txt' not in infilename:
		continue
	outfilename = infilename.replace('results.txt','tidy.results.txt')
	infile = open(indir + infilename,'r')
	outfile = open(indir + outfilename,'w')
	for line in infile:
		if 'ENSG' not in line:
			continue
		if 'rs' in line or 'mcgi' in line or 'ncgi' in line:
			continue
		liner = line.strip()
		print >> outfile, liner

	infile.close()
	outfile.close()
