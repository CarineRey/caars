#!/bin/bash

set -euo pipefail +o nounset

export HASH=`git rev-parse --short HEAD`
export BRANCH=`git branch | grep \* | cut -d ' ' -f2`

if [ $BRANCH = "master" ]
then
    ### CAARS_MASTER

    IMAGE_NAME=caars
    DOCKERFILE_DIR=caars

    export ENV_TAG="master_20200421"

    export TAG="master_$HASH"

else
    ### CAARS_DEV
    IMAGE_NAME=caars_dev
    DOCKERFILE_DIR=caars_dev

    export ENV_TAG="dev_"$BRANCH"_20200421"
    #export ENV_TAG="generax_20200331"

    export TAG="dev_$BRANCH_$HASH"
fi

REPO=carinerey/$IMAGE_NAME:$TAG

echo "## Build docker: $REPO ##"

docker build --no-cache --build-arg BRANCH=$BRANCH --build-arg ENV_TAG=$ENV_TAG -t $REPO $DOCKERFILE_DIR

push_flag=$1

if [[ $push_flag == "push_yes" ]]
then
    echo "## Push docker ##"
    docker push $REPO:$TAG
fi
