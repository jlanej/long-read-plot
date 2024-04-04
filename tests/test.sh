#!/bin/bash


# causes a script to immediately exit when it encounters an error.
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

bam=/long-read-plot/examples/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.bam
region=chr17:10958130-11017414

/app/NGS-TL/ngsTL.sh \
--cramFile "$cramFile" \
--craiFile "$craiFile" \
--rootOutput "$rootOutput" \
--referenceGenome "$referenceGenome" \
--gcBedFile "$gcBedFile" \
--regionsSearch "$regionsSearch"


diff "$rootOutput".ltl.estimate.txt.gz $(dirname $0)/reference.output.NA12878.ltl.estimate.txt.gz
status=$?
exit $status
