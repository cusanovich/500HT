rawmat = read.table("/mnt/lustre/home/cusanovich/oldhome_nb/Hutterite_Heritability/Idcoefs/addSNP.coef.3671",header=T)
orderlist = read.table("/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.rename.fam")[,2]

ordermat = matrix(NA,length(orderlist),length(orderlist))
for(i in 1:dim(rawmat)[1]){
	inds1 = match(rawmat[i,1],orderlist)
	inds2 = grep(rawmat[i,2],orderlist)
	ordermat[inds1,inds1] = rawmat[i,3]
	ordermat[inds2,inds1] = rawmat[i,3]
}

write.table(ordermat,"/mnt/lustre/home/cusanovich/500HT/addSNP.1415.ordered.txt",row.name=F,col.names=F,quote=F)
