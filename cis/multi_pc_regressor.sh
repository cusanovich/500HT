#!/bin/bash
for i in {17..20}
do
   nohup python /mnt/lustre/home/cusanovich/500HT/Scripts/combiner.py $i &
done