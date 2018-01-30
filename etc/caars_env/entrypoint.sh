#!/bin/bash

# Add local user
# Either use the LOCAL_USER_ID if passed in at runtime or
# fallback


USER_ID=${LOCAL_USER_ID:-9001}
Xvfb :1 -screen 0 1024x768x16 &
export DISPLAY=:1



echo "Working directory:" `pwd`
echo "(You can use -w docker option to change it)"

if [ -n "$LOCAL_USER_ID" ] 
then
USER_ID=${LOCAL_USER_ID:-9001}
echo "Starting with UID : $USER_ID"
echo "To be root, type \" sudo su - \""
useradd --shell /bin/bash -u $USER_ID -o -c "" -g sudo -m user_docker
echo "user_docker ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
export HOME=/home/user_docker
exec /usr/sbin/gosu user_docker "$@"
else
echo "Starting with UID : $USER_ID"
exec "$@"
fi
