 perl -lane 'print $1 if /(miR-\d+-*\d*\w*)/ig' RA-miRNA.txt | sort -u > RA_108_miRNA.txt
 tr '[:upper:]' '[:lower:]' < RA_108_miRNA.txt > RA_miRNA.txt
 sort -u RA_miRNA.txt > RA_miRNA.uni.txt
 
 perl -lane 'print $1 if /(miR-\d+)/ig' RA-miRNA.txt | sort -u > RA_108_miRNA.txt
 tr '[:upper:]' '[:lower:]' < RA_108_miRNA.txt > RA_miRNA.txt
 sort -u RA_miRNA.txt > RA_miRNA.uni.txt
