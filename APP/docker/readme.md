run docker run hello-world it gives me following error: WARNING: Error loading config file:/home/user/.docker/config.json - stat /home/user/.docker/config.json: permission denied. But then it continues and gives a success message:Hello from Docker.
This message shows that your installation appears to be working correctly.
```
sudo chown "$USER":"$USER" /home/"$USER"/.docker -R
sudo chmod g+rwx "/home/$USER/.docker" -R
sudo chown $(whoami):docker /home/$(whoami)/.docker/config.json
docker run hello-world
```
