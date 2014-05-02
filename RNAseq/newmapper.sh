##Script written by Bryce to map reads using BWA
##Adapted by DC 3/21/12, 3/24/12

REFGENOME="/mnt/lustre/data/share/HumanGenome/hg19/allhg19_norandom.fasta.gz"
REFJUNC="/mnt/lustre/home/bmvdgeijn/RNAseq/Mappers/all_junctions.50.ens.eedb"
#DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell1"
#DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell2"
#DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell3"
#DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell4"
#DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell5"
#DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell6"
#DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell7"
#DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell8"
#DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell9"
#DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell10"
DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell11"
#OUTDIR="/mnt/lustre/data/users/cusanovich/RNAseq_scripts/BWA_OE_Files/FlowCell1/"
#OUTDIR="/mnt/lustre/data/users/cusanovich/RNAseq_scripts/BWA_OE_Files/FlowCell2/"
#OUTDIR="/mnt/lustre/data/users/cusanovich/RNAseq_scripts/BWA_OE_Files/FlowCell3/"
#OUTDIR="/mnt/lustre/data/users/cusanovich/RNAseq_scripts/BWA_OE_Files/FlowCell4/"
#OUTDIR="/mnt/lustre/data/users/cusanovich/RNAseq_scripts/BWA_OE_Files/FlowCell5/"
#OUTDIR="/mnt/lustre/data/users/cusanovich/RNAseq_scripts/BWA_OE_Files/FlowCell6/"
#OUTDIR="/mnt/lustre/data/users/cusanovich/RNAseq_scripts/BWA_OE_Files/FlowCell7/"
#OUTDIR="/mnt/lustre/data/users/cusanovich/RNAseq_scripts/BWA_OE_Files/FlowCell8/"
#OUTDIR="/mnt/lustre/data/users/cusanovich/RNAseq_scripts/BWA_OE_Files/FlowCell9/"
#OUTDIR="/mnt/lustre/data/users/cusanovich/RNAseq_scripts/BWA_OE_Files/FlowCell10/"
OUTDIR="/mnt/lustre/data/users/cusanovich/RNAseq_scripts/BWA_OE_Files/FlowCell11/"
COUNTER="/mnt/lustre/data/users/cusanovich/RNAseq_scripts/counters.txt"

x=1
for NAME in `ls $DIR` # 106651 # 
	do
#	if [ $x -eq 16 ]||[ $x -gt 16 ]
#		then
#		break
#		fi
	echo $NAME
	for lane in `ls ${DIR}/${NAME}/*.sequence.txt.gz | sed 's/.saved//g' | uniq`
		do
		LANE_NAME=`echo $lane | sed 's/.sequence.txt.gz//g'` # | sed 's/\//./g'`
#		echo $lane
		echo $LANE_NAME
   
		if [ ! -f ${LANE_NAME}.quality.sort.bam ]
			then
#			echo Yipes!
			echo $x
#			echo $lane
			echo "\
			/mnt/lustre/home/cusanovich/Programs/bwa/bwa aln -n 2 -t 3  ${REFGENOME} ${lane} > ${LANE_NAME}.Ref.sai; \
			/mnt/lustre/home/cusanovich/Programs/bwa/bwa samse -n 1 ${REFGENOME} ${LANE_NAME}.Ref.sai ${lane} > ${LANE_NAME}.Ref.sam; \
			rm ${LANE_NAME}.Ref.sai; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -q 10 -b ${LANE_NAME}.Ref.sam > ${LANE_NAME}.quality.bam ; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -f 4 -b ${LANE_NAME}.Ref.sam > ${LANE_NAME}.unmapped.bam ; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools sort ${LANE_NAME}.quality.bam ${LANE_NAME}.quality.sort; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -c ${LANE_NAME}.Ref.sam > ${LANE_NAME}.count.txt; \
			rm ${LANE_NAME}.quality.bam; \
			rm ${LANE_NAME}.Ref.sam; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -c ${LANE_NAME}.quality.sort.bam >> ${LANE_NAME}.count.txt; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -c ${LANE_NAME}.unmapped.bam >> ${LANE_NAME}.count.txt; \
			/mnt/lustre/home/cusanovich/Programs/bwa/bwa aln -n 2 -t 3  ${REFJUNC} -b0 ${LANE_NAME}.unmapped.bam > ${LANE_NAME}.junction.Ref.sai; \
			/mnt/lustre/home/cusanovich/Programs/bwa/bwa samse -n 1 ${REFJUNC} ${LANE_NAME}.junction.Ref.sai ${LANE_NAME}.unmapped.bam > ${LANE_NAME}.junction.Ref.sam; \
			rm ${LANE_NAME}.junction.Ref.sai; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -q 10 -b ${LANE_NAME}.junction.Ref.sam > ${LANE_NAME}.junction.quality.bam ; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -f 4 -b ${LANE_NAME}.junction.Ref.sam > ${LANE_NAME}.junction.unmapped.bam ; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools sort ${LANE_NAME}.junction.quality.bam ${LANE_NAME}.junction.quality.sort; \
			rm ${LANE_NAME}.junction.quality.bam; \
			rm ${LANE_NAME}.junction.Ref.sam; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -c ${LANE_NAME}.junction.quality.sort.bam >> ${LANE_NAME}.count.txt; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -c ${LANE_NAME}.junction.unmapped.bam >> ${LANE_NAME}.count.txt; \
			/mnt/lustre/home/cusanovich/Programs/bwa/bwa aln -n 2 -t 3  ${REFGENOME} ${LANE_NAME}.saved.sequence.txt.gz > ${LANE_NAME}.saved.Ref.sai; \
			/mnt/lustre/home/cusanovich/Programs/bwa/bwa samse -n 1 ${REFGENOME} ${LANE_NAME}.saved.Ref.sai ${LANE_NAME}.saved.sequence.txt.gz > ${LANE_NAME}.saved.Ref.sam; \
			rm ${LANE_NAME}.saved.Ref.sai; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -q 10 -b ${LANE_NAME}.saved.Ref.sam > ${LANE_NAME}.saved.quality.bam ; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -f 4 -b ${LANE_NAME}.saved.Ref.sam > ${LANE_NAME}.saved.unmapped.bam ; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools sort ${LANE_NAME}.saved.quality.bam ${LANE_NAME}.saved.quality.sort; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -c ${LANE_NAME}.saved.Ref.sam >> ${LANE_NAME}.count.txt; \
			rm ${LANE_NAME}.saved.quality.bam; \
			rm ${LANE_NAME}.saved.Ref.sam; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -c ${LANE_NAME}.saved.quality.sort.bam >> ${LANE_NAME}.count.txt; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -c ${LANE_NAME}.saved.unmapped.bam >> ${LANE_NAME}.count.txt; \
			/mnt/lustre/home/cusanovich/Programs/bwa/bwa aln -n 2 -t 3  ${REFJUNC} -b0 ${LANE_NAME}.saved.unmapped.bam > ${LANE_NAME}.saved.junction.Ref.sai; \
			/mnt/lustre/home/cusanovich/Programs/bwa/bwa samse -n 1 ${REFJUNC} ${LANE_NAME}.saved.junction.Ref.sai ${LANE_NAME}.saved.unmapped.bam > ${LANE_NAME}.saved.junction.Ref.sam; \
			rm ${LANE_NAME}.saved.junction.Ref.sai; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -q 10 -b ${LANE_NAME}.saved.junction.Ref.sam > ${LANE_NAME}.saved.junction.quality.bam ; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -f 4 -b ${LANE_NAME}.saved.junction.Ref.sam > ${LANE_NAME}.saved.junction.unmapped.bam ; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools sort ${LANE_NAME}.saved.junction.quality.bam ${LANE_NAME}.saved.junction.quality.sort; \
			rm ${LANE_NAME}.saved.junction.quality.bam; \
			rm ${LANE_NAME}.saved.junction.Ref.sam; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -c ${LANE_NAME}.saved.junction.quality.sort.bam >> ${LANE_NAME}.count.txt; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -c ${LANE_NAME}.saved.junction.unmapped.bam >> ${LANE_NAME}.count.txt; \
			paste ${COUNTER} ${LANE_NAME}.count.txt > ${LANE_NAME}.counts.txt; \
			rm ${LANE_NAME}.count.txt" | qsub -l h_vmem=32g -l bigio=0 -o ${OUTDIR} -e ${OUTDIR} -wd `pwd` -N "map.${NAME}"
			x=$(( $x + 1 ))
		fi
		done
	done
#			/mnt/lustre/home/cusanovich/Programs/samtools/samtools merge ${LANE_NAME}.total.quality.sort.bam ${LANE_NAME}.quality.sort.bam ${LANE_NAME}.junction.quality.sort.bam ${LANE_NAME}.saved.quality.sort.bam ${LANE_NAME}.saved.junction.quality.sort.bam"; \
