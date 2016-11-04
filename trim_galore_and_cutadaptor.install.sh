
# install pip
wget https://bootstrap.pypa.io/get-pip.py

pip install --user --upgrade cutadapt
~/.local/bin/cutadapt --help
cp ~/.local/bin/cutadapt ~/bin/
cp ~/.local/bin/cutadapt /media/Home_Raid1/shg047/software/trim_galore_zip

git clone https://github.com/marcelm/cutadapt.git
bismark -q --phred33-quals -n 1 -l 40 --non_directional ~/NAS3/
bismark -q --phred33-quals -n 1 -l 30 ~/NAS3/db/hg19/ -1 T84_R1.fastq -2 T84_R2.fastq
bismark -q --phred33-quals -n 1 -l 30 ~/NAS3/db/hg19/ -1 T84_R1_val_1.fq -2 T84_R2_val_2.fq

T84_R1_val_1.fq

