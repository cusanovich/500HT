exprs = read.table("/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.newcovcor.bimbam.gz")
#exprs = read.table("~/home/500HT/qqnorm.500ht.bimbam")

htpca = prcomp(exprs,scale.=TRUE)
xhtpca = htpca$x
write.table(xhtpca,"/mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.newcovcor.pcs",
            col.names=F,row.names=F,quote=F,sep="\t")
