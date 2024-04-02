library(optparse)
source("bamRUtils.R")
# options for the bam input file, image output file, and the region to plot in ucsc format
option_list = list(
  make_option(
    c("-i", "--bam"),
    type = "character",
    default = "./examples/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.bam",
    help = "bam file to plot",
    metavar = "file"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "./examples/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.png",
    help = "output image file",
    metavar = "file"
  ),
  make_option(
    c("-r", "--region"),
    type = "character",
    default = "chr17:10958130-11017414",
    help = "region to plot in ucsc format",
    metavar = "region"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
bamFile = opt$bam
# parse the region to plot to Granges
region = strsplit(opt$region, ":|-")[[1]]
bamAll = parseAlignments(bamFile, region)
# write bamAll to a file named opt$output with a txt extension and the output extension removed
write.table(
  apply(bamAll, 2, as.character),
  file = paste0(tools::file_path_sans_ext(opt$output), ".txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
