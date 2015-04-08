#!/usr/bin/env python
import gzip
infile = gzip.open('/mnt/lustre/home/cusanovich/snp138.txt.gz','rb')
outfile = gzip.open('/mnt/lustre/home/cusanovich/snp138.cropped.txt.gz','wb')
keepcols = [1,2,3,4,12,16,17]
x = 0
for line in infile:
	liner = line.strip().split('\t')
	if liner[17] != '1':
		continue
	if liner[18] != '':
		continue
	print >> outfile, "\t".join([liner[x] for x in keepcols])
	x += 1

print str(x) + " records retained."
infile.close()
outfile.close()
