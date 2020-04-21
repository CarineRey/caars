#!/bin/bash

set -euo pipefail +o nounset

export HASH=`git rev-parse --short HEAD`
export BRANCH=`git branch | grep \* | cut -d ' ' -f2`
export DATE=`date +'%Y%m%d'`


IMAGE_NAME=caars_env
DOCKERFILE_DIR=caars_env

echo BRANCH : $BRANCH

if [ $BRANCH = "master" ] 
then
    TAG="master_$DATE"
else
    TAG="dev_"$BRANCH"_"$DATE
fi

REPO=carinerey/$IMAGE_NAME:$TAG

echo "## Build docker: $REPO ##"

cp -r ../utils $DOCKERFILE_DIR
docker build -t $REPO $DOCKERFILE_DIR
rm -r $DOCKERFILE_DIR/utils

push_flag=$1

if [[ $push_flag == "push_yes" ]]
then
    echo "## Push docker ##"
    docker push $REPO
fi
