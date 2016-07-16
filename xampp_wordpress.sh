# install xampp
sudo chmod 777 -R xampp-linux-x64-5.6.23-0-installer.run
sudo ./xampp-linux-x64-5.6.23-0-installer.run
# stop and restart
sudo /etc/init.d/apache2 stop
sudo /etc/init.d/mysql stop
sudo /etc/init.d/proftpd stop
sudo /opt/lampp/lampp start
sudo service vsftpd stop
sudo /opt/lampp/lampp stop
# start xampp
sudo /opt/lampp/lampp start
# check
http://132.239.189.199/

# install wordpress(bitnami)
sudo chmod 777 -R bitnami-wordpress-3.9.1-1-module-linux-x64-installer.run
sudo ./bitnami-wordpress-3.9.1-1-module-linux-x64-installer.run
# set blog name (Kun Zhang's lab, email, passwd)

# ifconfig to find IP
http://132.239.189.199/wordpress/

