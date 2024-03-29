#!/bin/bash


singularity run --bind "$(pwd):$(pwd)" "docker://ghcr.io/jlanej/long-read-plot:main"