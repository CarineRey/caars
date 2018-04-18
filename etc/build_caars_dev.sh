#!/bin/bash

set -euo pipefail +o nounset

### CAARS_DEV
IMAGE_NAME=caars_dev
DOCKERFILE_DIR=caars_dev
REPO=carinerey/$IMAGE_NAME
export BRANCH_DEV=`git branch | grep \* | cut -d ' ' -f2`
TAG=$BRANCH_DEV
echo BRANCH_DEV=$BRANCH_DEV
docker build  --no-cache --build-arg BRANCH_DEV="$BRANCH_DEV" -t $REPO:$TAG $DOCKERFILE_DIR
push_flag=$1

if [[ $push_flag == "push_yes" ]]
then
    docker push $REPO:$TAG
fi
