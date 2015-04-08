exprs = read.table("/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.newcovcor.bimbam.gz")
exprnames = read.table("/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.newcovcor.findivs.txt")
pcs = read.table("/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.newcovcor.pcs")
idcoefs = read.table("/mnt/lustre/home/cusanovich/500HT/addSNP.coef.3671.square")
matcher = read.table("/mnt/lustre/home/cusanovich/500HT/3chip/hutt.3chip.fam")
matcherids = matcher$V2
matcherind = match(matcherids,exprnames[,1])
exprs.o = exprs[matcherind,]
pcs.o = pcs[matcherind,]
idcoefs.o = idcoefs[matcherind,matcherind]
write.table(pcs.o,"/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.newcovcor.3chip_order.pcs",row.names=F,col.names=F,quote=F,sep="\t")
write.table(exprs.o,"/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.newcovcor.3chip_order.bimbam",row.names=F,col.names=F,quote=F,sep="\t")
write.table(idcoefs.o,"/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.newcovcor.3chip_order.square.txt",row.names=F,col.names=F,quote=F,sep="\t")

