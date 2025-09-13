#!/bin/bash


PYTHON=/venv-sbml2sim/bin/python3

docker run -it \
    -e DISPLAY=$DISPLAY --net=host \
    --name sbml2sim2 sbml2sim2 \
    $PYTHON main.py 

docker rm -f sbml2sim2 2>/dev/null || true