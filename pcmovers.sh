#!/bin/bash
for i in $(seq 0 1 40)
do
   cp /mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.gccor.covcor.3chip_order.pc$i /mnt/lustre/home/cusanovich/500HT/Exprs/qqnorm.500ht.gccor.covcor.ordered.pc$i
done
