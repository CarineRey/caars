#!/bin/bash

set -eo

### CAARS
IMAGE_NAME=caars
DOCKERFILE_DIR=caars
REPO=carinerey/$IMAGE_NAME
docker build --no-cache -t $REPO $DOCKERFILE_DIR
docker push carinerey/$IMAGE_NAME
