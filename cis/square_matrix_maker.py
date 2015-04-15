print "Building coefficient dictionary..."
rawmat = open("/mnt/lustre/home/cusanovich/oldhome_nb/Hutterite_Heritability/Idcoefs/addSNP.coef.3671","r")
covdic = {}
for line in rawmat:
    liner = line.strip().split()
    covdic[(liner[0],liner[1])] = liner[2]

rawmat.close()
print "Making new matrix..."
orderfile = open("/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.500ht.fam","r")
orderlist = [x.strip().split()[1] for x in orderfile.readlines()]
orderfile.close()
outfile = open("/mnt/lustre/home/cusanovich/500HT/addSNP.500ht.ordered.square.test.txt","w")
for inda in orderlist:
    currline = ['NA']*len(orderlist)
    currcount = 0
    for indb in orderlist:
        ind1 = inda
        ind2 = indb
        if int(inda) > int(indb):
            ind1 = indb
            ind2 = inda
        currline[currcount] = covdic[(ind1,ind2)]
        currcount += 1
    print >> outfile, "\t".join(currline)

outfile.close()

