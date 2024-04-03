#!/bin/bash

bam=$HOME/git/long-read-plot/examples/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.bam
output=/scratch.global/lanej/1000G/long_read/plots/NA19240_output/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.png
region=chr17:10958130-11017414

singularity run --containall \
--pwd "/long-read-plot/" \
--bind "$HOME/git/long-read-plot/examples" \
--bind "/scratch.global/lanej/1000G/long_read/plots/NA19240_output" \
"docker://ghcr.io/jlanej/long-read-plot:main" \
Rscript longReadPlot.R --bam $bam --output $output --region $region --debug