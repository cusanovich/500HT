library(scales)
eqtls = read.table("500HT/eQTLs/master.PC62.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.chosen.txt")
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
h = hist(distance[sigind],breaks=50,plot=F)
h1 = hist(distance[-sigind],breaks=50,plot=F)

pdf("~/eqtl_distance.pdf")
hist(distance[-sigind],breaks=50,col=alpha("indianred",0.5),freq=F,ylim=c(0,5e-06),
     main="Distance from TSS for cis-eQTLs",xlab="Distance from TSS",las=1)
hist(distance[sigind],breaks=50,col=alpha("dodgerblue2",0.5),add=T,freq=F)
legend("topright",c("Nonsig. cis-eQTLs","Sig. cis-eQTLs"),
       fill=c(alpha("indianred",0.5),alpha("dodgerblue2",0.5)))
plot(density(distance[sigind]),col="dodgerblue2",xlab="Distance from TSS",lwd=2,
     main="Distance from TSS for cis-eQTLs",xaxp=c(-1000000,1000000,8),las=1)
lines(density(distance[-sigind]),col="indianred",lwd=2)
abline(v=0,lty="dashed",lwd=2)
legend("topright",c("Nonsig. cis-eQTLs","Sig. cis-eQTLs"),
       fill=c("indianred","dodgerblue2"))
plot(x=h$mids, y=h$density, type="l",las=1,xaxt="n",lwd=3,col="dodgerblue2",ylab="Normalized Density",xlab="Distance from TSS")
axis(1,at=seq(-1000000,1000000,250000))
lines(x=h1$mids, y=h1$density,lwd=3,col="indianred")
abline(v=0,lty="dashed",lwd=2)
legend("topright",c("Nonsig. cis-eQTLs","Sig. cis-eQTLs"),
       fill=c("indianred","dodgerblue2"))
dev.off()
