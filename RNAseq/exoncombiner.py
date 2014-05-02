"""Take in an exon count file and write out a gene count file."""
import sys
import operator

# This doesn't seem to work
#def gene_sorter(dicto):
#    """Sort a dictionary of value lists by value 0, then value 1
#
#    Source:
#    http://writeonly.wordpress.com/2008/08/30/sorting-dictionaries-by-value-in-python-improved/
#    &
#    http://wiki.python.org/moin/HowTo/Sorting
#
#    """
#    sort_start = sorted(dicto.iteritems(), key = operator.itemgetter(1))
#    return sorted(sort_start, key = operator.itemgetter(0))
#    return sort_start

f = open(sys.argv[1],'r')
genes = {}
for line in f:
    if line == "":
        continue
    else:
        line = line.strip().split('\t')
#        print line
#        break
        gene = line[3]
        if gene not in genes.keys():
            genes[gene] = [line[0],int(line[1]),int(line[2]),float(line[4]),int(line[6])]
            #genes[gene] = [line[0],int(line[1]),int(line[2]),int(line[4]),int(line[6])]
        else:
            genes[gene][1] = min(int(line[1]),genes[gene][1])
            genes[gene][2] = max(int(line[2]),genes[gene][2])
            #genes[gene][3] += int(line[4])
            genes[gene][3] += float(line[4])
            genes[gene][4] += int(line[6])
#        print liner
#        break
f.close()

#genes_sorted = gene_sorter(genes)

outfile = open(sys.argv[2],'w')
print >> outfile, "Chr\tStart\tEnd\tENSGID\tCount\tGeneLength"
#for key in genes_sorted:
for key in sorted(genes.iterkeys()):
    print >> outfile, genes[key][0] + '\t' + str(genes[key][1]) + '\t' + str(genes[key][2]) + '\t' + key + '\t' + str(genes[key][3]) + '\t' + str(genes[key][4])
#    print >> outfile, genes_sorted[key]
#    print >> outfile, key
outfile.close()
