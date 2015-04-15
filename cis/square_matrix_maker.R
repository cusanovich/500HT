rawmat = read.table("/mnt/lustre/home/cusanovich/oldhome_nb/Hutterite_Heritability/Idcoefs/addSNP.coef.3671",header=T)
orderlist = read.table("/mnt/lustre/home/cusanovich/500HT/Imputed1415/hutt.imputed.rename.fam")[,2]

orderorderlist = sort(orderlist)
newmat = matrix(NA,length(orderlist),length(orderlist))
for(i in 1:length(orderlist)){
  for(j in i:length(orderlist)){
    row_index = rawmat[, 1] == orderorderlist[i] & rawmat[, 2] == orderorderlist[j]
    newmat[i, j] = rawmat[row_index, 3]
    newmat[j, i] = rawmat[row_index, 3]  
  }
}

unorderind = match(orderlist,orderorderlist)
ordermat = newmat[unorderind,unorderind]
write.table(ordermat,"/mnt/lustre/home/cusanovich/500HT/addSNP.1415.ordered.txt",row.name=F,col.names=F,quote=F)
