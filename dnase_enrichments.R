print("Starting the engines...")
dnase = read.table("/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.dnase.bed", header=T)
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
pdf("/mnt/lustre/home/cusanovich/500HT/Imputed1415/tissues2.pdf")
for(i in 1:length(gwass)){
	gwas = read.table(paste0("/mnt/lustre/home/cusanovich/500HT/MinalSherlock/",gwass[i]),header=T)
	currpheno = strsplit(gwass[i],"[.]")[[1]][1]
	print(currpheno)
	currind = match(gwas[,1],dnase[,4])
	currdnase = dnase[currind,5:dim(dnase)[2]]
	noners = which(is.na(currdnase[,1]))
	currdnase[noners,] = rep(0,dim(currdnase)[2])
	baselines = colSums(currdnase>0)/dim(currdnase)[1]
	currthresholds = quantile(gwas[,2],thresholds)
	currenrichments = matrix(NA,dim(currdnase)[2],length(currthresholds))
	print("Calculating enrichments...")
	for(j in 1:dim(currdnase)[2]){
		for(k in 1:length(currthresholds)) {
			currenrichments[j,k] = length(which(gwas[,2] < currthresholds[k] & currdnase[,j] > 0))/length(which(gwas[,2] < currthresholds[k]))/baselines[j]
		}
	}
	currenrichments = cbind(currenrichments,rep(1,times=dim(currdnase)[2]))
	colnames(currenrichments) = c(currthresholds,1)
	rownames(currenrichments) = colnames(currdnase)
	write.table(currenrichments,paste0("/mnt/lustre/home/cusanovich/500HT/MinalSherlock/",currpheno,".enrichmenttable.txt"),quote=F,sep="\t")
	print("Plotting enrichments...")
	plot(-log10(as.numeric(colnames(currenrichments))),currenrichments[1,],type="l",ylim=c(0,max(currenrichments)),main=currpheno,col=cellcolor[1])
	for(i in 2:dim(currenrichments)[1]){
		lines(-log10(as.numeric(colnames(currenrichments))),currenrichments[i,],col=cellcolor[i])
	}
}
dev.off()
