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
exprs = matrix_reader('/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.newcovcor.3chip_order.bimbam',sep=" ",dtype='|S10')
#mastercols = matrix_reader('/mnt/lustre/home/cusanovich/500HT/hutt.3chip.500kb.mastercols.txt',dtype='|S15')

generecs = open('/mnt/lustre/home/cusanovich/500HT/ensemblCAGETSS_RNAseq_sorted.bed','r')
genedic = {}
for line in generecs:
	genedic[line.strip().split()[3]] = line.strip().split()[0]

mastergenes = open('/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.newcovcor.genenames.txt','r')
genes = mastergenes.readlines()
genes = [x.strip() for x in genes]


print "Building dictionaries..."
masterdic = {}
exprcoldic = {}
chrmdic = {}
#for i in range(mastercols.shape[0]):
for col,gene in enumerate(genes):
	try:
		#masterdic[mastercols[i,0]].append(mastercols[i,1])
		masterdic[gene].append(col)
	except KeyError:
		masterdic[gene] = [col]
		exprcoldic[gene] = col
		chrmdic[gene] = genedic[gene]

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
	numpy.savetxt('/mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.gccor.newcovcor.ordered.chr' + str(chrm) + '.bimbam',currexpr,delimiter=" ",fmt='%s')
	genelist = open('/mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.gccor.newcovcor.ordered.chr' + str(chrm) + '.genes','w')
	print >> genelist, "\n".join(genes)
	genelist.close()
