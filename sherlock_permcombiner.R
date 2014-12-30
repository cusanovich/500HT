permfile = read.table("/mnt/lustre/home/cusanovich/500HT/Imputed1415/perms/output/perm_1.assoc.txt",header=T)

permmat = matrix(NA,dim(permfile)[1],51)
permmat[,1] = as.character(permfile[,2])
permmat[,2] = permfile[,9]

for(i in 2:50){
	print(i)
	currfile = read.table(paste0("/mnt/lustre/home/cusanovich/500HT/Imputed1415/perms/output/perm_",i,".assoc.txt"),header=T)
	ind = match(permmat[,1],currfile[,2])
	permmat[,(i+1)] = currfile[ind,9]
	permmat = permmat[!is.na(permmat[,(i+1)]),]
}

print(dim(permmat))
write.table(permmat,"/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.permfile",row.names=F,col.names=F,quote=F,sep="\t")
