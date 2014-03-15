transers = read.table("~/500HT/trans_maximizer.txt")

pdf("~/500HT/trans_maximizer.pdf")
plot(c(0:59),log10(transers$V2),pch=20,ylim=c(2.5,3.5),xlab = "No. PCs Removed",
     ylab = "Log10(No. Significant Associations Between TFs and Targets)",
     type="b",col="dodgerblue2")
points(c(0:59),log10(transers$V4),type="b",col="indianred",pch=20)
legend("topright",c("Knockdown Targets","Random Targets"),
       fill=c("dodgerblue2","indianred"))
plot(c(0:59),transers$V1,pch=20,xlab = "No. PCs Removed",
     ylab = "Median Correlation between TFs and Targets",type="b",col="dodgerblue2")
points(c(0:59),transers$V3,pch=20,type="b",col="indianred")
legend("topright",c("Knockdown Targets","Random Targets"),
       fill=c("dodgerblue2","indianred"))
dev.off()