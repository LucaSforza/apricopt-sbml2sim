#!/bin/bash

if [ -n "$1" ]; then
    DOCKERFILE="$1"
else
    DOCKERFILE="Dockerfile.base"
fi

docker build -f $DOCKERFILE --platform linux/amd64 --pull --rm -t sbml2sim2:base .