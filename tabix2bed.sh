#!/bin/bash
for i in {1..22}
do
   zcat /mnt/lustre/home/cusanovich/500HT/Imputed1415/ByChr/hutt.imputed.chr$i.txt.gz | cut -f1-4 | gzip > /mnt/lustre/home/cusanovich/500HT/Imputed1415/ByChr/hutt.imputed.chr$i.bed.gz
done
