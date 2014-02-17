#! /usr/bin/Rscript --vanilla
library(EDASeq)
args <- commandArgs(TRUE)

gctable = read.table("/mnt/lustre/home/cusanovich/500HT/exon.gccontent.txt",header=T)
expressedgenes = read.table("/mnt/lustre/home/cusanovich/500HT/genenames.500ht.txt")
exoncounts = read.table(args[1])
outfile = gsub("exoncounts","gccor.exoncounts",args[1])

expressed.ind = c()
for(i in 1:dim(expressedgenes)[1]){
  expressed.ind = c(expressed.ind,grep(expressedgenes[i,1],gctable[,1]))
}

expressing = sort(unique(expressed.ind))
gctable.expr = gctable[expressing,]
exoncounts.expr = exoncounts[expressing,]
cleanedup = withinLaneNormalization(exoncounts.expr,gctable.expr,which="full",round=F)

write.table(cleanedup,outfile,row.names=F,col.names=F,quote=F,sep="\t")