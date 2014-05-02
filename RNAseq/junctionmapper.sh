##Script written by Bryce to map reads using BWA
##Adapted by DC 2/10/12 to remap the unmapped reads to an exon junction db

REF1a="/mnt/lustre/home/bmvdgeijn/RNAseq/Mappers/all_junctions.50.ens.eedb"
DIR="/mnt/lustre/data/users/cusanovich/500HTRNAseq/FlowCell2"

x=1
for NAME in `ls $DIR` #106651
	do
	if [ $x -eq 32 ]||[ $x -gt 32 ]
		then
		break
		fi
	echo $NAME
	for lane in `ls ${DIR}/${NAME}/*.unmapped.bam | uniq`
		do
		LANE_NAME=`echo $lane | sed 's/.unmapped.bam//g'`
		echo $LANE_NAME
   
		if [ ! -f ${LANE_NAME}.junction.quality.sort.bam ]&&[[ ! $LANE_NAME = *.junction ]]
			then
#			echo Yipes!
			echo $x
#			echo $lane
#			echo $LANE_NAME
			echo "\
			/mnt/lustre/home/cusanovich/Programs/bwa/bwa aln -n 2 -t 3  ${REF1a} -b0 ${LANE_NAME}.unmapped.bam > ${LANE_NAME}.junction.Ref.sai; \
			/mnt/lustre/home/cusanovich/Programs/bwa/bwa samse -n 1 ${REF1a} ${LANE_NAME}.junction.Ref.sai ${LANE_NAME}.unmapped.bam > ${LANE_NAME}.junction.Ref.sam; \
			rm ${LANE_NAME}.junction.Ref.sai; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -q 10 -b ${LANE_NAME}.junction.Ref.sam > ${LANE_NAME}.junction.quality.bam ; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -S -h -f 4 -b ${LANE_NAME}.junction.Ref.sam > ${LANE_NAME}.junction.unmapped.bam ; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools sort ${LANE_NAME}.junction.quality.bam ${LANE_NAME}.junction.quality.sort; \
			rm ${LANE_NAME}.junction.quality.bam; \
			rm ${LANE_NAME}.junction.Ref.sam; \
			/mnt/lustre/home/cusanovich/Programs/samtools/samtools view -c ${LANE_NAME}.junction.quality.sort.bam >> ${LANE_NAME}.count.txt" | qsub -l h_vmem=8g -l bigio=0 -wd `pwd` -N "juncmap"
			x=$(( $x + 1 ))
		fi
		done
	done