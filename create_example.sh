#!/bin/bash


PYTHON=/venv-sbml2sim/bin/python3

input_file=$1
output_dir=$2

docker run -it \
    -e DISPLAY=$DISPLAY --net=host \
    --name sbml2sim2 sbml2sim2:base $PYTHON sbml2sim/generate_experiment.py $input_file $output_dir

docker cp sbml2sim2:/app/$output_dir .

docker rm -f sbml2sim2 2>/dev/null || true