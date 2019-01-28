# search 1092 RA genes + IL7 and IL7R
# and then find level 1 and level2 network

grep IL7 IL-7_RA.tsv > output.txt
grep CD4 IL-7_RA.tsv >> output.txt
awk '{print $1}' output.txt > output2.txt
awk '{print $2}' output.txt >> output2.txt
sort -u output2.txt > output3.txt
