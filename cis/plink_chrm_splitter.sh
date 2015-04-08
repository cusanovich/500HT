#!/bin/bash
for i in {1..22}
do
	plink --bfile /mnt/lustre/home/cusanovich/500HT/3chip/hutt.3chip --chr $i --make-bed --out /mnt/lustre/home/cusanovich/500HT/3chip/ByChr/chr$i.hutt.3chip
done
