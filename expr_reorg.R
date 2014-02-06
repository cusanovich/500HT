exprs = read.table("/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.bimbam")
exprnames = read.table("/mnt/lustre/home/cusanovich/500HT/findivs.500ht.txt")
pcs = read.table("/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.pcs")
matcher = read.table("/mnt/lustre/home/cusanovich/500HT/3chip/hutt.3chip.fam")
matcherids = matcher$V2
matcherind = match(matcherids,exprnames[,1])
exprs.o = exprs[matcherind,]
pcs.o = pcs[matcherind,]
write.table(pcs.o,"/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.3chip_order.pcs",row.names=F,col.names=F,quote=F,sep="\t")
write.table(exprs.o,"/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.3chip_order.bimbam",row.names=F,col.names=F,quote=F,sep="\t")

