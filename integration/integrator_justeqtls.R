gwass = list.files("/mnt/lustre/home/cusanovich/500HT/MinalSherlock", pattern="sherlock")
eqtls = read.table("/mnt/lustre/home/cusanovich/500HT/eQTLs/master.PC62.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.chosen.txt")
converter = read.table("/mnt/lustre/home/cusanovich/500HT/sherlock_converter.txt")
eqtlq = p.adjust(eqtls[,4],method="BH")

justeqtls = eqtls[which(eqtlq < 0.05),]

perming = function(currgwas,lister){
	permgwas = currgwas[sample(c(1:dim(currgwas)[1])),]
	permmeds = sapply(lister,function(x) median(permgwas[1:x,2]))
	return(permmeds)
}

lister = c(1:1000)
pdf("/mnt/lustre/home/cusanovich/justeqtl_pval_enrichments.pdf")
for(i in 1:length(gwass)){
	currpheno = strsplit(gwass[i],"[.]")[[1]][1]
	print(currpheno)
	print("Cleaning data...")
	gwas = read.table(paste0("/mnt/lustre/home/cusanovich/500HT/MinalSherlock/",gwass[i]))
	exprs = read.table(paste0("/mnt/lustre/home/cusanovich/500HT/dege/output/",converter[match(gwass[i],converter[,1]),2]),header=T)
	commongenes = match(intersect(exprs[,2],justeqtls[,1]),justeqtls[,1])
	commonsnps = match(intersect(gwas[,1],justeqtls[,2]),justeqtls[,2])
	curreqtls = justeqtls[intersect(commongenes,commonsnps),]
	currgwasind = match(curreqtls[,2],gwas[,1])
	currgwas = gwas[currgwasind,]
	currexprind = match(curreqtls[,1],exprs[,2])
	currexpr = exprs[currexprind,c(2,13)]
	currrank = order(currexpr[,2])
	currgwas = currgwas[currrank,]
	currexpr = currexpr[currrank,]
	grandmed = median(currgwas[,2])
	currmeds = sapply(lister,function(x) median(currgwas[1:x,2]))
	print("Permuting data...")
	currperms = replicate(1000,perming(currgwas,lister))
	ups = apply(currperms,1,function(x) quantile(x,0.95))
	downs = apply(currperms,1,function(x) quantile(x,0.05))
	print("Plotting data...")
	plot(currmeds,ylim=c(0,1),type="l",main=currpheno,col="indianred",lwd=2)
	abline(h=grandmed,lty="dashed",lwd=2)
	lines(c(1:1000),ups,lwd=2,col="dodgerblue2")
	lines(c(1:1000),downs,lwd=2,col="dodgerblue2")
}
dev.off()
