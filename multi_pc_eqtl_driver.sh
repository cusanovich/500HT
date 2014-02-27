#!/bin/bash
for i in $(seq 0 1 40)
do
   nohup python /mnt/lustre/home/cusanovich/500HT/Scripts/eqtl_driver.py $i &
done
