
Step 1: install GSL
wget http://mirror.keystealth.org/gnu/gsl/gsl-latest.tar.gz
tar xzvf gsl-latest.tar.gz
cd /home/shg047/software/gsl-2.1
./configure --prefix=/media/Home_Raid1/shg047/software/gsl-2.2.1
make
make install
export CPATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/include
export LIBRARY_PATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/lib

Step2： install methpipe
wget http://smithlabresearch.org/downloads/methpipe-3.4.2.tar.bz2
tar xjvf methpipe-3.4.2.tar.bz2
make
﻿export PATH=$PATH:/media/Home_Raid1/shg047/software/methpipe-3.4.2/bin

Step 3: install rmap
wget http://smithlabresearch.org/downloads/rmap-2.1.tar.bz2
tar xjvf rmap-2.1.tar.bz2
cd rmap-2.1/
make
sudo make install
PATH=$PATH:/media/Home_Raid1/shg047/software/rmap-2.1/bin

Step 4. Install walt
git clone https://github.com/smithlabcode/walt.git
cd walt
make

Step 5: Download Genome Reference - hg19
cd /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz


Step 6: 
makedb -c /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome -o methpipe.hg19.dbindex



