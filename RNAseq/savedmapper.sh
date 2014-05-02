##Script written by Bryce to map reads using BWA
##Adapted by DC 3/21/12

REF1a=/mnt/lustre/data/share/HumanGenome/hg19/allhg19_norandom.fasta.gz
DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell2"

x=1
for NAME in `ls $DIR` #106651 #
	do
#	if [ $x -eq 16 ]||[ $x -gt 16 ]
#		then
#		break
#		fi
	echo $NAME
	for lane in `ls ${DIR}/${NAME}/*.saved.sequence.txt.gz | sed 's/.unmapped//g' | uniq`
		do
		LANE_NAME=`echo $lane | sed 's/.sequence.txt.gz//g'` # | sed 's/\//./g'`
		echo $LANE_NAME
		LANER=`echo $lane | sed 's/.saved.sequence.txt.gz//g'` # | sed 's/\//./g'`
   
		if [ ! -f ${LANE_NAME}.quality.sort.bam ]
			then
#			echo Yipes!
			echo $x
			echo $lane
			echo "\
			/mnt/lustre/home/cusanovich/Programs/bwa/bwa aln -n 2 -t 3  ${REF1a} ${lane} > ${LANE_NAME}.Ref.sai; \
			/mnt/lustre/home/cusanovich/Programs/bwa/bwa samse -n 1 ${REF1a} ${LANE_NAME}.Ref.sai ${lane} > ${LANE_NAME}.Ref.sam; \
			rm ${LANE_NAME}.Ref.sai; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -q 10 -b ${LANE_NAME}.Ref.sam > ${LANE_NAME}.quality.bam ; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -f 4 -b ${LANE_NAME}.Ref.sam > ${LANE_NAME}.unmapped.bam ; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools sort ${LANE_NAME}.quality.bam ${LANE_NAME}.quality.sort; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -c ${LANE_NAME}.Ref.sam >> ${LANER}.count.txt; \
			rm ${LANE_NAME}.quality.bam; \
			rm ${LANE_NAME}.Ref.sam; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -c ${LANE_NAME}.quality.sort.bam >> ${LANER}.count.txt; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -c ${LANE_NAME}.unmapped.bam >> ${LANER}.count.txt; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools merge ${LANER}.total.quality.sort.bam ${LANER}.quality.sort.bam ${LANER}.saved.quality.sort.bam" | qsub -l h_vmem=8g -l bigio=0 -wd `pwd` -N "savemap.${NAME}"
			x=$(( $x + 1 ))
		fi
		done
	done

		#	/mnt/lustre/home/cusanovich/Programs/bwa/bwa aln -n 2 -t 3  ${REF1a} -b0 ${LANE_NAME}.unmapped.bam > ${LANE_NAME}.junction.Ref.sai; \
		#	/mnt/lustre/home/cusanovich/Programs/bwa/bwa samse -n 1 ${REF1a} ${LANE_NAME}.junction.Ref.sai ${LANE_NAME}.unmapped.bam > ${LANE_NAME}.junction.Ref.sam; \
		#	rm ${LANE_NAME}.junction.Ref.sai; \
		#	/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -q 10 -b ${LANE_NAME}.junction.Ref.sam > ${LANE_NAME}.junction.quality.bam ; \
		#	/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -f 4 -b ${LANE_NAME}.junction.Ref.sam > ${LANE_NAME}.junction.unmapped.bam ; \
		#	/mnt/lustre/home/cusanovich/Programs/samtools/samtools sort ${LANE_NAME}.junction.quality.bam ${LANE_NAME}.junction.quality.sort; \
		#	rm ${LANE_NAME}.junction.quality.bam; \
		#	rm ${LANE_NAME}.junction.Ref.sam; \
		#	/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -c ${LANE_NAME}.junction.quality.sort.bam >> ${LANER}.count.txt; \