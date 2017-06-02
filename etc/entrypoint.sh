#!/bin/bash

# Add local user
# Either use the LOCAL_USER_ID if passed in at runtime or
# fallback

if [ -n "$SHARED_DIR" ]
then
W_DIR=$SHARED_DIR
fi

W_DIR=${W_DIR:-/data}

echo "Working directory : $W_DIR"
cd $W_DIR

if [ -n "$LOCAL_USER_ID" ] 
then
USER_ID=${LOCAL_USER_ID:-9001}
echo "Starting with UID : $USER_ID"
useradd --shell /bin/bash -u $USER_ID -o -c "" -g sudo -m user_amalgam
adduser user_amalgam sudo
echo "user_amalgam ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
export HOME=/home/user_amalgam
exec /usr/local/bin/gosu user_amalgam "$@"
else
echo "Starting with UID : root"
exec "$@"
fi
