#!/bin/bash

# Add local user
# Either use the LOCAL_USER_ID if passed in at runtime or
# fallback

USER_ID=${LOCAL_USER_ID:-9001}
W_DIR=${W_DIR:-/data}


echo "Starting with UID : $USER_ID"
cd $W_DIR
echo "Working directory : $W_DIR"
useradd --shell /bin/bash -u $USER_ID -o -c "" -m user_amalgam
export HOME=/home/user_amalgam

exec /usr/local/bin/gosu user_amalgam "$@"
