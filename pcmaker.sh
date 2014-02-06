#!/bin/bash
for i in {1..20}
do
   paste <(cat /mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.3chip_order.pc0) <(cut -f1-$i /mnt/lustre/home/cusanovich/500HT/qqnorm.500ht.3chip_order.pcs) > /mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.3chip_order.pc$i
done
