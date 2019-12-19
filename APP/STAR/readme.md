### STAR mapping pipeline with 2-pass for multiple samples
Shicheng Guo, 2019/12/18, San Francisco

STAR Installation-1
```
ubuntu.
$ sudo apt-get update
$ sudo apt-get install g++
$ sudo apt-get install make
Red Hat, CentOS, Fedora.
$ sudo yum update
$ sudo yum install make
$ sudo yum install gcc-c++
$ sudo yum install glibc-static
SUSE.
$ sudo zypper update
$ sudo zypper in gcc gcc-c++
```
STAR Installation-2
```
cd ~/hpc/tools/
wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz
tar -xzf 2.7.3a.tar.gz
cd STAR-2.7.3a
```
Indexing genome with annotations
```
cd  ~/hpc/db/
wget http://ftp.ensemblorg.ebi.ac.uk/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
gunzip Homo_sapiens.GRCh37.87.gtf.gz
STAR --runMode genomeGenerate --genomeDir ~/hpc/db/hg19/STAR/ --genomeFastaFiles ~/hpc/db/hg19/hg19.fa --sjdbGTFfile ~/hpc/db/hg19/Homo_sapiens.GRCh37.87.gtf --runThreadN 30 --sjdbOverhang 89
```


Reference Download Website(RDW):

GTF: http://useast.ensembl.org/info/data/ftp/index.html
