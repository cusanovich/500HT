exonfiles=`find /mnt/lustre/home/cusanovich/data/500HTRNAseq/ -name "*gccor.exoncounts.txt" -exec echo {} \;`
echo "Exon files collected..." 
for exonfile in $exonfiles
	do
		genefile=`echo $exonfile | sed 's/exoncounts/genecounts/g'`
		echo "python /mnt/lustre/home/cusanovich/data/RNAseq_scripts/exoncombiner.py $exonfile $genefile" | qsub -l h_vmem=2g -V -o ~/dump/ -e ~/dump/ -N "gene_count" 
	done
