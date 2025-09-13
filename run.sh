#!/bin/bash


PYTHON=/venv-sbml2sim/bin/python3

docker rm -f sbml2sim2 2>/dev/null || true

docker run -it \
    -e DISPLAY=$DISPLAY --net=host \
    --name sbml2sim2 sbml2sim2:base bash

docker rm -f sbml2sim2 2>/dev/null || true