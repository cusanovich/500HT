import gzip
infile = gzip.open('/mnt/lustre/home/cusanovich/snp138.cropped.txt.gz','rb')
outfile = open('/mnt/lustre/home/cusanovich/snp138.unique.txt','w')

baddic = {}
for line in infile:
	liner = line.strip().split('\t')
	badkey = ':'.join(liner)
	try:
		baddic[badkey].append(liner[0:3])
	except KeyError:
		print >> outfile, '\t'.join(liner)
		baddic[badkey] = [liner]

infile.close()
outfile.close()

badfile = outfile = open('/mnt/lustre/home/cusanovich/snp138.AlsoMulti.txt','w')
for snp in baddic.keys():
	for record in snp:
		print >> badfile, '\t'.join(record)

badfile.close()
