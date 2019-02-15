```
cd
wget https://github.com/ldc-developers/ldc/releases/download/v$ver/ldc2-1.7.0-linux-x86_64.tar.xz
tar xvJf ldc2-1.7.0-linux-x86_64.tar.xz
export PATH=$HOME/ldc2-1.7.0-linux-x86_64/bin:$PATH
export LIBRARY_PATH=$HOME/ldc2-1.7.0-linux-x86_64/lib
git clone --recursive https://github.com/biod/sambamba.git
cd sambamba
make
```
