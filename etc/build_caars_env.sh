#!/bin/bash

set -eo

IMAGE_NAME=caars_env
DOCKERFILE_DIR=caars_env
TAG="0.2.0-bistrodev"
REPO=carinerey/$IMAGE_NAME:$TAG
cp -r ../utils $DOCKERFILE_DIR
docker build -t $REPO $DOCKERFILE_DIR
rm -r $DOCKERFILE_DIR/utils

if [[ $1 == "push_yes" ]]
then
    docker push $REPO
fi
