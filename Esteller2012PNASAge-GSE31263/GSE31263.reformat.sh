#!/usr/bin/sh
for i in `ls *bed`
do
perl GSE31263.reformat.pl $i > $i.bedGraph
done
