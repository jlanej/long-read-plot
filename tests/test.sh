#!/bin/bash


# causes a script to immediately exit when it encounters an error.
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

bam=/long-read-plot/examples/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.bam
outputDir=$HOME/tmp/NA19240_output/
outputRoot=$outputDir/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414
region=chr17:10958130-11017414

Rscript /long-read-plot/longReadPlot.R --bam $bam --output $output.png --region $region --debug

#files created by the debug run above
filesCreated=( "$outputRoot".adjusted.tsv.gz "$outputRoot".png "$outputRoot".tsv.gz )

# files to compare are located in the image here
compareDir=/long-read-plot/tests/test_output
filesToCompare=($compareDir/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.adjusted.tsv.gz $compareDir/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.png $compareDir/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.tsv.gz)

# diff the files and exit if they are different
for i in ${!filesCreated[@]}; do
    diff ${filesCreated[$i]} ${filesToCompare[$i]}
    status=$?
    if [ $status -ne 0 ]; then
        exit $status
    fi
done
