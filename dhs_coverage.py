infile = open("/mnt/lustre/home/cusanovich/500HT/dhs_and_gaps_union.bed","r")
outfile1 = open("/mnt/lustre/home/cusanovich/500HT/dhs_union_sizes.txt","w")
outfile2 = open("/mnt/lustre/home/cusanovich/500HT/dhs_gap_union_sizes.txt","w")
chrm = ''
gaps = []
for line in infile:
	liner = line.strip().split()
	if liner[3] == '1':
		print >> outfile1, str(int(liner[2]) - int(liner[1]))
		if liner[0] != chrm:
			chrm = liner[0]
		continue
	if liner[3] == '0':
		if liner[0] == 'chrY' and liner[1] == '59033141':
			print >> outfile2, "\n".join(gaps)
			continue
		if liner[0] == chrm:
			gaps.append(str(int(liner[2]) - int(liner[1])))
		if liner[0] != chrm and len(gaps) > 0:
			print >> outfile2, "\n".join(gaps[0:(len(gaps)-1)])
			gaps = []
			continue
