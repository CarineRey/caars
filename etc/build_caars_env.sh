#!/bin/bash

set -euo pipefail +o nounset

IMAGE_NAME=caars_env
DOCKERFILE_DIR=caars_env
TAG="master"
REPO=carinerey/$IMAGE_NAME:$TAG
cp -r ../utils $DOCKERFILE_DIR
docker build -t $REPO $DOCKERFILE_DIR
rm -r $DOCKERFILE_DIR/utils

push_flag=$1

if [[ $push_flag == "push_yes" ]]
then
    docker push $REPO
fi
