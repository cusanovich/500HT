#!/bin/bash
for i in $(seq 90 10 181)
do
   nohup python /mnt/lustre/home/cusanovich/500HT/Scripts/eqtl_driver.py $i &
done
