print("Starting the engines...")
dnase = read.table("/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.dnase2.bed", header=T)
gwass = list.files("/mnt/lustre/home/cusanovich/500HT/MinalSherlock", pattern="sherlock")
thresholds = c(0.00025,0.0005,0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128,0.256,0.5)
coloring = c("indianred", "darkred", "red", "dodgerblue2",
             "aquamarine", "cadetblue1", "mediumseagreen", "darkgreen",
             "yellowgreen", "mediumorchid3", "violetred2", "darkmagenta",
             "orange3", "coral4", "orangered2", "gold", "yellow3")
celltype = read.table("/mnt/lustre/home/cusanovich/500HT/stam_cells.txt")
cellcolors = coloring[unlist(as.factor(celltype[,2]))]
colorind = match(colnames(dnase)[5:length(colnames(dnase))],celltype[,1])
cellcolor = cellcolors[colorind]
badcells = which(is.na(cellcolor))
cellcolored = cellcolor[-badcells]
dnased = dnase[,-(badcells+4)]
celltype = NA
dnase = NA
pdf("/mnt/lustre/home/cusanovich/tissue_enrichments.pdf")
for(i in 1:length(gwass)){
	gwas = read.table(paste0("/mnt/lustre/home/cusanovich/500HT/MinalSherlock/",gwass[i]))
	currpheno = strsplit(gwass[i],"[.]")[[1]][1]
	print(currpheno)
	currind = match(gwas[,1],dnased[,4])
	currdnase = dnased[currind,5:dim(dnased)[2]]
	noners = which(is.na(currdnase[,1]))
	currdnase[noners,] = rep(0,dim(currdnase)[2])
	baselines = colSums(currdnase>0)/dim(currdnase)[1]
	currthresholds = quantile(gwas[,2],thresholds)
	currenrichments = matrix(NA,dim(currdnase)[2],length(currthresholds))
	currps = rep(1,dim(currdnase)[2])
	print("Calculating enrichments...")
	for(j in 1:dim(currdnase)[2]){
		for(k in 1:length(currthresholds)){
			currenrichments[j,k] = length(which(gwas[,2] < currthresholds[k] & currdnase[,j] > 0))/length(which(gwas[,2] < currthresholds[k]))/baselines[j]
			if(k == 1){
				olaps = length(which(gwas[,2] < currthresholds[k] & currdnase[,j] > 0))
				dnasers = length(which(currdnase[,j] > 0))
				nondnasers = dim(currdnase)[1] - dnasers
				chancers = length(which(gwas[,2] < currthresholds[k]))
				currps[j] = phyper((olaps - 1),dnasers,nondnasers,chancers, lower.tail=FALSE)
			}
		}
	}
	currqs = p.adjust(currps,method="BH")
	currenrichments = cbind(currenrichments,rep(1,times=dim(currdnase)[2]))
	colnames(currenrichments) = c(currthresholds,1)
	rownames(currenrichments) = colnames(currdnase)
	#write.table(currenrichments,paste0("/mnt/lustre/home/cusanovich/500HT/MinalSherlock/",currpheno,".enrichmenttable.txt"),quote=F,sep="\t")
	print("Plotting enrichments...")
	plot(-log10(as.numeric(colnames(currenrichments))),currenrichments[1,],type="l",ylim=c(0,max(currenrichments)),main=currpheno,col="gray")
	for(l in 2:dim(currenrichments)[1]){
		lines(-log10(as.numeric(colnames(currenrichments))),currenrichments[l,],col="gray")
	}
	for(m in 1:dim(currenrichments)[1]){
		if(currqs[m] < 0.2){
                	lines(-log10(as.numeric(colnames(currenrichments))),currenrichments[m,],col=cellcolored[m])
		}
        }
	currenrichments = cbind(currqs,currenrichments)
	colnames(currenrichments)[1] = "Qvals"
        currenrichments = cbind(currps,currenrichments)
        colnames(currenrichments)[1] = "Pvals"
	write.table(currenrichments,paste0("/mnt/lustre/home/cusanovich/500HT/MinalSherlock/",currpheno,".enrichmenttable.txt"),quote=F,sep="\t")
}
dev.off()
