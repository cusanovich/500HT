setwd("/mnt/lustre/home/cusanovich/500HT/dege/")
phcvts = list.files("./",pattern=".ph-cvt")
forder = read.table("../findivs.500ht.ordered.txt")
pdf("/mnt/lustre/home/cusanovich/covariate_plots.pdf")
par(mfrow=c(2,2))
for(i in 1:length(phcvts)){
#for(i in 1:4){
	pheno = read.table(phcvts[i],header=T)
        phenoname = gsub("[.]ph-cvt","",phcvts[i])
	print(phenoname)
	porder = match(forder[,1],pheno[,1])
	ages = grep("age",colnames(pheno))
	inbreeds = grep("inbreed",colnames(pheno))
	hts = grep("ht",colnames(pheno))
	HTs = grep("HT",colnames(pheno))
	quantcons = c(ages,inbreeds,hts,HTs)
	factcons = 4:dim(pheno)[2]
	factcons = setdiff(factcons,quantcons)
	if(length(factcons) > 0){
		for(j in 1:length(factcons)){
			print(factcons[j])
			boxplot(pheno[porder,2] ~ as.factor(pheno[porder,factcons[j]]),main=phenoname,xlab=colnames(pheno)[factcons[j]])
		}
	}
	if(length(quantcons) > 0){
		for(k in 1:length(quantcons)){
			print(quantcons[k])
			plot(pheno[porder,quantcons[k]],pheno[porder,2],pch=20,xlab=colnames(pheno)[quantcons[k]],ylab=phenoname)
		}
	}
}
dev.off()
