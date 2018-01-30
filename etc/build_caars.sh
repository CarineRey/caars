#!/bin/bash

set -euo pipefail +o nounset

### CAARS
IMAGE_NAME=caars
DOCKERFILE_DIR=caars

export HASH=`git rev-parse --short HEAD`
TAG=$HASH
REPO=carinerey/$IMAGE_NAME:$TAG

docker build --no-cache --build-arg BRANCH_DEV=$BRANCH -t $REPO $DOCKERFILE_DIR
push_flag=$1

if [[ $push_flag == "push_yes" ]]
then
    docker push $REPO
fi
