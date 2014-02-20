#! /usr/bin/Rscript --vanilla
library(aroma.light,lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.0/")
library(EDASeq,lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.0/")
args = commandArgs(TRUE)
#args = c("~/data/500HTRNAseq/FlowCell1/106272/lane_3.index_8.exoncounts.txt")
print("Loading tables...")
gctable = read.table("/mnt/lustre/home/cusanovich/500HT/exon.gccontent.txt",header=T)
expressedgenes = read.table("/mnt/lustre/home/cusanovich/500HT/genenames.500ht.txt")
exoncounts = read.table(args[1])
outfile = gsub("exoncounts","gccor.exoncounts",args[1])
print("GC correcting...")
expressing = gctable[,1] %in% expressedgenes[,1]
gctable.expr = gctable[expressing,]
exoncounts.expr = exoncounts[expressing,]
cleanedup = withinLaneNormalization(as.matrix(exoncounts.expr[,5]),as.vector(gctable.expr[,3]),which="full",round=F)
exoncounts.cleanedup = exoncounts.expr
exoncounts.cleanedup$V5 = cleanedup
print("Writing results...")
write.table(exoncounts.cleanedup,outfile,row.names=F,col.names=F,quote=F,sep="\t")
