setwd("/mnt/lustre/home/cusanovich/500HT/dege/")
phcvts = list.files("./",pattern="_phenotype.txt")
for(i in 1:length(phcvts)){
#for(i in 1:4){
	pheno = read.table(phcvts[i])
	phenoname = gsub("_phenotype.txt","",phcvts[i])
	print(phenoname)
	system(paste0("~/Programs/gemma0.94 -p ", phenoname, "_phenotype.txt -k ../addSNP.500ht.ordered.square.txt -c ", phenoname, "_covariates.txt -lmm 5 -o ", phenoname,"_regress"))
	regre = read.table(paste0("output/",phenoname,"_regress.residE.txt"))
	notna = which(!is.na(pheno[,1]))
	regreind = order(regre[,1])
	regremid = regreind[11:(length(regreind)-10)]
	regrehi = regreind[1:10]
	regrelo = regreind[(length(regreind) - 9):length(regreind)]
	#regremid = regreind[c((round(length(regreind)/2,0)-10):(round(length(regreind)/2,0)+10))]
	#pheno[notna[regremid],1] = NA
	pheno[notna[regrelo],1] = NA
	pheno[notna[regrehi],1] = NA
	write.table(pheno,paste0(phenoname,"_phenotype_extreme.txt"),col.names=F,row.names=F,quote=F)
	system(paste0("~/Programs/gemma0.94 -g expression_gemma_format.txt -p ", phenoname, "_phenotype_extreme.txt -k ../addSNP.500ht.ordered.square.txt -c ", phenoname, "_covariates.txt -notsnp -lmm 4 -o ", phenoname,"_extreme"))
	tester = read.table(paste0("output/",phenoname,"_extreme.assoc.txt"),header=T)
	print(min(p.adjust(tester[,13],method="BH")))
	print(length(which(p.adjust(tester[,13],method="BH") < 0.2)))
}
