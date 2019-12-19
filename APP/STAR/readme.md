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
Indexing genome with annotations (pay attention to chrosome name: 1,2,3 or chr1,chr2,chr3)
```
# https://www.gencodegenes.org/
cd  ~/hpc/db/hg19
wget http://ftp.ensemblorg.ebi.ac.uk/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/gencode.v29lift37.annotation.gtf.gz
wget wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/gencode.v32lift37.annotation.gtf.gz

gunzip Homo_sapiens.GRCh37.87.gtf.gz
STAR --runMode genomeGenerate --genomeDir ~/hpc/db/hg19/STAR/ --genomeFastaFiles ~/hpc/db/hg19/hg19.fa --sjdbGTFfile ~/hpc/db/hg19/gencode.v32lift37.annotation.gtf.gz --runThreadN 30 --sjdbOverhang 89
```


Reference Download Website(RDW):

GTF: http://useast.ensembl.org/info/data/ftp/index.html and [wget link](http://ftp.ensemblorg.ebi.ac.uk/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz)

FASTA:http://useast.ensembl.org/info/data/ftp/index.html and [wget link](http://ftp.ensemblorg.ebi.ac.uk/pub/grch37/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.alt.fa.gz)
