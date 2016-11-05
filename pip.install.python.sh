mkdir pip
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py --prefix="./pip"
pip install --user --upgrade cutadapt
