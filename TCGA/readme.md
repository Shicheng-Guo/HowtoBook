## The Cancer Genome Atlas (TCGA) Research Network

#### Full dataset for Causal Network Analysis

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
