#!/bin/bash

bam=/app/long-read-plot/examples/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.bam

singularity run  --pwd "/app/long-read-plot/" "docker://ghcr.io/jlanej/long-read-plot:main" Rscript longReadPlot.R

