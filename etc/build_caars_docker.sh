#!/bin/bash

set -euo pipefail +o nounset

export HASH=`git rev-parse --short HEAD`
export BRANCH=`git branch | grep \* | cut -d ' ' -f2`
push_flag=$1
branch_flag=$2

if [[ ! -z $branch_flag ]]
then
    #For travis
    echo "Replace BRANCH for travis"
    echo "from: $BRANCH"
    echo "to: $branch_flag"
    BRANCH=$branch_flag
fi

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

    #export ENV_TAG="dev_"$BRANCH"_20200421"
    export ENV_TAG="generax_20200331"

    export TAG="dev_$BRANCH_$HASH"
fi

REPO=carinerey/$IMAGE_NAME:$TAG

if [[ $push_flag == "only_pull" ]]
then
    echo "## Pull docker $REPO ##"
    docker pull $REPO
    exit 0
fi

if [[ $push_flag == "only_pull_env" ]]
then
    ENV_DEPO=carinerey/caars_env:$ENV_TAG
    echo "## Pull docker $ENV_DEPO ##"
    docker pull $ENV_DEPO
    exit 0
fi

echo "## Build docker: $REPO ##"
docker build --no-cache --build-arg BRANCH=$BRANCH --build-arg ENV_TAG=$ENV_TAG -t $REPO $DOCKERFILE_DIR

if [[ $push_flag == "push_yes" ]]
then
    echo "## Push docker ##"
    docker push $REPO
fi
