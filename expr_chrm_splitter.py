import subprocess
import numpy

def matrix_reader(matrix_file,sep="\t",dtype='|S20'):
	linecounter = subprocess.Popen('wc -l ' + matrix_file, shell=True, stdout=subprocess.PIPE)
	linecount = int(linecounter.communicate()[0].strip().split()[0])
	columncounter = subprocess.Popen('awk -F"' + sep + '" \'{print NF;exit}\' ' + matrix_file, shell=True, stdout=subprocess.PIPE)
	columncount = int(columncounter.communicate()[0].strip().split()[0])
	raws = numpy.zeros((linecount,columncount),dtype=dtype)
	rawin = open(matrix_file,'r')
	for i,line in enumerate(rawin):
		raws[i,:] = line.strip().split()
	rawin.close()
	return raws

print "Loading master files..."
exprs = matrix_reader('/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.3chip_order.bimbam',sep=" ",dtype='|S10')
mastercols = matrix_reader('/mnt/lustre/home/cusanovich/500HT/hutt.3chip.500kb.mastercols.txt',dtype='|S15')

print "Building dictionaries..."
masterdic = {}
exprcoldic = {}
chrmdic = {}
for i in range(mastercols.shape[0]):
	try:
		masterdic[mastercols[i,0]].append(mastercols[i,1])
	except KeyError:
		masterdic[mastercols[i,0]] = [mastercols[i,1]]
		exprcoldic[mastercols[i,0]] = mastercols[i,2]
		chrmdic[mastercols[i,0]] = mastercols[i,4]

"Print writing out files..."
for chrm in range(1,23):
	print chrm
	indices = []
	genes = []
	for gene in exprcoldic.keys():
		if chrmdic[gene] == 'chr' + str(chrm):
			indices.append(exprcoldic[gene])
			genes.append(gene)
	currexpr = exprs[:,indices]
	numpy.savetxt('/mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.3chip_order.chr' + str(chrm) + '.bimbam',currexpr,delimiter=" ",fmt='%s')
	genelist = open('/mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.3chip_order.chr' + str(chrm) + '.genes','w')
	print >> genelist, "\n".join(genes)
	genelist.close()
