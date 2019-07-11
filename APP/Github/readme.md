1. install and push local repository to Github 
```
sudo apt-get install git
git config --global user.name "Shicheng-Guo"
git config --global user.email "Shicheng.Guo@hotmail.com"
git init Mytest
Initialized empty Git repository in /home/akshay/Mytest/.git/
cd Mytest
gedit README
git add .
git commit -m "some_message"
git remote add origin https://github.com/user_name/Mytest.git
git push 
```

2. How to avoid input username and passwd for `git push`
```
cd ~
ssh-keygen -t rsa 
```
3. Copy the key to your github [setting page: https://github.com/settings/keys](https://github.com/settings/keys)

```
ssh-keygen -t rsa -b 4096 -C "Shicheng.Guo@hotmail.com"
ssh-add -K ~/.ssh/id_rsa
ssh-add -A
vim ~/.ssh/config

ls ~/.ssh/
sudo apt-get install xclip
xclip -sel clip < ~/.ssh/id_rsa.pub

eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_rsa

git remote set-url origin git@github.com:Shicheng-Guo/LungMethBiomarker.git
```


