#!/bin/bash
baselayer='/mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.gccor.covcor.3chip_order.pc0'
#baselayer='/mnt/lustre/home/cusanovich/500HT/knowncovariates.500ht.3chip_order.txt'
for i in $(seq 1 1 40)
do
   paste <(cat $baselayer) <(cut -f1-$i /mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.gccor.covcor.3chip_order.pcs) > /mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.gccor.covcor.3chip_order.pc$i
done
