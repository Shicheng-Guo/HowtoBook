Protocol To Downlad TCGA Data From GDC in Windows and to network directory
```
Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\luad

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
