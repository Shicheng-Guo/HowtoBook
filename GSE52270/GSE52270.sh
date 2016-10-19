for i in `ls GSM*`
do
perl GSE52270.pl $i > $i.bedgraph &
done



