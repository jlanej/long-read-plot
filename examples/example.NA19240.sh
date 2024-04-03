#!/bin/bash

#This example bam is packaged within image
bam=/long-read-plot/examples/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.bam
region=chr17:10958130-11017414
outputDir=$HOME/tmp/NA19240_output/
output=$outputDir/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.png

singularity run --containall \
--pwd "/long-read-plot/" \
--bind "$outputDir/" \
"docker://ghcr.io/jlanej/long-read-plot:main" \
Rscript longReadPlot.R --bam $bam --output $output --region $region --debug

