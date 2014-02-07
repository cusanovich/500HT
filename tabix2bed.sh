#!/bin/bash
for i in {1..22}
do
   zcat /mnt/lustre/home/cusanovich/500HT/3chip/ByChr/hutt.3chip.chr$i.txt.gz | cut -f1-4 | gzip > /mnt/lustre/home/cusanovich/500HT/3chip/ByChr/hutt.3chip.chr$i.bed.gz
done
