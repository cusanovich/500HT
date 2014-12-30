sherlockdir = "/mnt/lustre/home/cusanovich/500HT/MinalSherlock/"
sherlockers = list.files(sherlockdir,pattern=".tidy.results.txt")
for(i in 1:length(sherlockers)){
	currsher = read.table(paste0(sherlockdir,sherlockers[i]))
	adjps = p.adjust(currsher[,3],method="BH")
	print(sherlockers[i])
	print(length(which(adjps < 0.15)))
}
