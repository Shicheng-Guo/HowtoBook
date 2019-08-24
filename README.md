# HowtoBook

How to install R-3.3.0 in RedHat EL6.0
#### install anaconda3 for python 3.xxx
#### install anaconda2 for python 2.xxx
#### Be careful, some software can only be run in python 2 or python 3, for example, BSSeerker2 only works in python 2 rather than python 

3. Update (Python vs R)

#### install pip by anaconda3 (pip in anaconda3 maybe not updated, you can update with )
pip install --upgrade pip

#### install Cython by anaconda3 
pip install Cython --install-option="--no-cython-compile"

#### pysam (pysam have been installed in python 3.5, now want to install python 2.7)
#### install pysam by pip
pip uninstall pysam
pip install pysam
#### install pysam by conda
conda config --add channels r
conda config --add channels bioconda
conda install pysam
