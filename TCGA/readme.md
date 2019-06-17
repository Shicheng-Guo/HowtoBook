## The Cancer Genome Atlas (TCGA) Research Network

#### Full dataset for Causal Network Analysis

Marshfield Clinic don't support Linux gdc-client. I download gdc-client.exe and run it in Windows System. 

LUAD: Lung Adenocarcinoma (TCGA-LUAD) data collection 
```
cd /gpfs/home/guosa/hpc/tcga/xiong/luad

mkdir mh450
mkdir fpkm-uq
mkdir miRNA
mkdir mCNS
mkdir SNP

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\luad\mh450
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\luad\fpkm-uq
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\luad\miRNA
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\luad\mCNS
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\luad\clinic
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\luad\SNP
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt
```


HNSC: Head-Neck Squamous Cell Carcinoma (TCGA-HNSC) data collection 
```
cd /gpfs/home/guosa/hpc/tcga/xiong/hnsc 

mkdir mh450
mkdir fpkm-uq
mkdir miRNA
mkdir mCNS
mkdir SNP

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\mh450
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\fpkm-uq
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\miRNA
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\mCNS
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\clinic
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\SNP
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt
```
Sample Size for each cancer for DNA methylation 450K dataset
```
	1	11
KIRC	324	160
BRCA	791	96
MESO	87	87
THCA	507	56
HNSC	528	50
PRAD	502	50
LIHC	377	50
UCEC	438	46
KIRP	275	45
LUSC	370	42
COAD	313	38
LUAD	473	32
BLCA	418	21
ESCA	185	16
PAAD	184	10
OV.H	10	10
CHOL	36	9
READ	98	7
SARC	261	4
PCPG	179	3
CESC	307	3
GBM.	140	2
SKCM	105	2
STAD	395	2
THYM	124	2
LGG.	516	0
TGCT	150	0
UVM.	80	0
ACC.	80	0
KICH	66	0
UCS.	57	0
DLBC	48	0
```
BRCA+LUNG
```
cd /gpfs/home/guosa/hpc/methylation/BRCA/450K/beta
```

