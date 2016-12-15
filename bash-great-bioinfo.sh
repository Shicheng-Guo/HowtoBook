#How to run GREAT in the batch mode to 100-1000 bed files at the same time.
#Our genome-miner server: http://132.239.25.238/shg047/
#bed files: /opt/lampp/htdocs/shg047/NAS3/Alice/WGBS/permutation/
#check the apache status,if it is okay: http://132.239.25.238/hello.pl
#Then run the following code in genome-miner:
#cd /opt/lampp/htdocs/shg047/NAS3/Alice/WGBS/permutation/hg19

for i in `ls *bed`
do
genome="hg19"
wget -O $i.results.tsv "http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?outputType=batch&requestSpecies=$genome&requestName=Example+Data&requestSender=Client+A&requestURL=http%3A%2F%2F132.239.25.238%2Fshg047%2FNAS3%2FAlice%2FWGBS%2Fpermutation%2Fhg19%2F$i"
done
cd /opt/lampp/htdocs/shg047/NAS3/Alice/WGBS/permutation/mm9
for i in `ls *bed`
do
genome="mm9"
echo wget -O $i.results.tsv "http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?outputType=batch&requestSpecies=$genome&requestName=Example+Data&requestSender=Client+A&requestURL=http%3A%2F%2F132.239.25.238%2Fshg047%2FNAS3%2FAlice%2FWGBS%2Fpermutation%2Fmm9%2F$i"
done
