for i in {2..50}
do
	#echo $i
	echo "python /mnt/lustre/home/cusanovich/500HT/Scripts/sherlock_reorder_genos.py $i" | qsub -l h_vmem=2g -l bigio=0 -o ~/dump -e ~/dump -wd `pwd` -N "genoperm.${i}" 
done
