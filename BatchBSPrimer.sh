#/usr/bin/sh
bedtools getfasta -fi ~/oasis/db/hg19_lambda.fa -bed venn.bed -fo esca.fa.bsp.fa
perl -p -i -e 's/CG/X/ig' esca.fa.bsp.fa
perl -p -i -e 's/C/T/ig' esca.fa.bsp.fa
perl -p -i -e 's/X/<CG>/ig' esca.fa.bsp.fa
perl -p -i -e 's/Thr/chr/ig' esca.fa.bsp.fa
# send esca.fa.bsp.fa to http://batchprimer3.bioinformatics.ucdavis.edu/cgi-bin/batchprimer3/batchprimer3.cgi 

bedtools getfasta -fi ~/oasis/db/hg19_lambda.fa -bed venn.t100.bed -fo esca.fa.bsp.fa
perl -p -i -e 's/CG/X/ig' esca.fa.bsp.fa
perl -p -i -e 's/C/T/ig' esca.fa.bsp.fa
perl -p -i -e 's/X/<CG>/ig' esca.fa.bsp.fa
perl -p -i -e 's/Thr/chr/ig' esca.fa.bsp.fa


