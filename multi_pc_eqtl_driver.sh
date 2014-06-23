#!/bin/bash
for i in $(seq 61 1 100)
do
   nohup python /mnt/lustre/home/cusanovich/500HT/Scripts/eqtl_driver.py $i &
done
