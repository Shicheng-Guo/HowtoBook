1. install and push local repository to Github 
```
sudo apt-get install git
git config --global user.name "user_name"
git config --global user.email "email_id"
git init Mytest
Initialized empty Git repository in /home/akshay/Mytest/.git/
cd Mytest
gedit README
git add .
git commit -m "some_message"
git remote add origin https://github.com/user_name/Mytest.git
git push origin master
```
