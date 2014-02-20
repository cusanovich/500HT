#library(qvalue)
library(plyr)
library(stringr)

all.pvals <- list.files(path = "/mnt/lustre/home/cusanovich/500HT/eQTLs/",pattern=".3chip.150kb.bonferroni.gemma.chosen.txt",full.names=T)
pvals <- llply(all.pvals, read.table)

pval.pcs.names = str_split(all.pvals,"[.]")
pval.pcs = c()
for(i in 1:length(pval.pcs.names)){
  pval.pcs[i] = str_split(pval.pcs.names[[i]][2],"PC")[[1]][2]
}
names(pvals) = pval.pcs
inder = order(as.numeric(pval.pcs))
pvals.o = pvals[inder]
pdf("~/500HT/Scripts/hists.pdf")
qcounts = c()
for(i in 1:length(pvals.o)){
  #qs = qvalue(pvals[[i]]$V4)$qvalues
  qs = p.adjust(pvals.o[[i]]$V4,method="BH")
  qcount = sum(qs < 0.05)
  qcounts[i] = qcount
  print(qcount)
  hist(pvals.o[[i]]$V4,xlab="P-value",main=paste0(names(pvals.o)[i]," PCs Removed"))
}
dev.off()
pdf("~/500HT/Scripts/eqtl.calls.pdf")
plot(as.numeric(names(pvals.o)),qcounts,xlab="Number of PCs Removed",ylab="No. of eQTLs",type="b",pch=20)
dev.off()
