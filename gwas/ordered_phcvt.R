inphcvt = read.table("/mnt/lustre/home/cusanovich/500HT/dege/lnIgeCheckSNP.ph-cvt",header=T)
outphcvt = "/mnt/lustre/home/cusanovich/500HT/dege/lnIgeCheckSNP_full.ph-cvt"
outbimbam = "/mnt/lustre/home/cusanovich/500HT/dege/lnIgeCheckSNP_full.pheno.txt"
outcov = "/mnt/lustre/home/cusanovich/500HT/dege/lnIgeCheckSNP_full.covariates.txt"
orderlist = read.table("/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.rename.fam")[,2]

updatemat = matrix(NA,length(orderlist),dim(inphcvt)[2])
ind = match(inphcvt[,1],orderlist)
for(i in 1:dim(inphcvt)[2]){
	updatemat[ind,i] = inphcvt[,i]
}
colnames(updatemat) = colnames(inphcvt)
updatemat[,1] = orderlist

write.table(updatemat,outphcvt,row.names=F,quote=F,sep="\t")
write.table(as.matrix(updatemat[,2]),outbimbam,row.names=F,col.names=F,quote=F)
write.table(updatemat[,3:dim(updatemat)[2]],outcov,row.names=F,col.names=F,quote=F,sep="\t")
