#!/usr/bin/bash
for i in `ls GSE17972_HUMtg5lib.qmap*.txt`
do
perl qmap2bedgraph.pl $i > $i.bedgraph &
done

