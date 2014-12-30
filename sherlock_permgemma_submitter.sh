for i in {2..50}
do
	cd /mnt/lustre/home/cusanovich/500HT/Imputed1415/perms/;
	echo "~/Programs/gemma0.94 -g /mnt/lustre/home/cusanovich/500HT/Imputed1415/perms/perm.${i}.genos.bimbam -p /mnt/lustre/home/cusanovich/500HT/dege/lnIgeCheckSNP_full.pheno.txt -k /mnt/lustre/home/cusanovich/500HT/addSNP.1415.ordered.txt -c /mnt/lustre/home/cusanovich/500HT/dege/lnIgeCheckSNP_full.covariates.txt -lmm 2 -o perm_${i}" | qsub -l h_vmem=4g -l bigio=0 -o ~/dump -e ~/dump -wd `pwd` -N "gemmaperm.${i}"
done
