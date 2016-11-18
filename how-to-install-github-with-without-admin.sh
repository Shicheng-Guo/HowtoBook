# How to install git with adminstration 
sudo apt-get update
sudo apt-get install git
git config --global user.name "Shicheng-Guo"
git config --global user.email "shicheng.guo@hotmail.com"
git config --list
git commit --amend --reset-author

# How to install git without admin
cd
mkdir git
cd git
wget https://www.kernel.org/pub/software/scm/git/git-2.10.2.tar.xz
tar xf xzvf git-2.10.2.tar.xz
cd git-2.10.2
./configure --prefix="$HOME/git"
make install
git config --list
