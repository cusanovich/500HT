DIR="/mnt/lustre/home/cusanovich/500HT/ByChr/"
for file in `ls ${DIR} | grep 150kbbonferroni | grep chosen`
do
#echo $DIR$file
filenew=`echo $file | sed 's/150kbbonferroni/150kb.bonferroni/g'`
#echo $filenew
mv $DIR$file $DIR$filenew
done