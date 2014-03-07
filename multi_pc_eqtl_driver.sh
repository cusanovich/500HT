#!/bin/bash
for i in $(seq 41 1 60)
do
   nohup python /mnt/lustre/home/cusanovich/500HT/Scripts/eqtl_driver.py $i &
done
