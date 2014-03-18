library(scales)
eqtls = read.table("500HT/eQTLs/master.PC42.imputed.150kb.bonferroni.gccor.covcor.regressPCs.gemma.chosen.txt")
tss = read.table("500HT/ensemblCAGETSS_RNAseq_sorted.bed")
snps = read.table("500HT/Imputed1415/hutt.imputed.500ht.bim")
geneind = match(eqtls$V1,tss$V4)
snpind = match(eqtls$V2,snps$V2)
distance = tss[geneind,3] - snps[snpind,4]
multi = gsub("-",-1,tss$V6)
multi = gsub("[+]",1,multi)
multi = as.numeric(multi)
distance = distance*multi[geneind]
sigind = which(p.adjust(eqtls$V4,method="BH") < 0.05)

pdf("eqtl_distance.pdf")
hist(distance[-sigind],breaks=25,col=alpha("indianred",0.5),freq=F,ylim=c(0,1e-05),
     main="Distance from TSS for cis-eQTLs",xlab="Distance from TSS")
hist(distance[sigind],breaks=25,col=alpha("dodgerblue2",0.5),add=T,freq=F)
legend("topright",c("Nonsig. cis-eQTLs","Sig. cis-eQTLs"),
       fill=c(alpha("indianred",0.5),alpha("dodgerblue2",0.5)))
plot(density(distance[sigind]),col="dodgerblue2",xlab="Distance from TSS",lwd=2,
     main="Distance from TSS for cis-eQTLs",xaxp=c(-150000,150000,6))
lines(density(distance[-sigind]),col="indianred",lwd=2)
abline(v=0,lty="dashed",lwd=2)
legend("topright",c("Nonsig. cis-eQTLs","Sig. cis-eQTLs"),
       fill=c("indianred","dodgerblue2"))
dev.off()
