library(optparse)
source("bamRUtils.R")
# options for the bam input file, image output file, and the region to plot in ucsc format
# also an option for debug
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
    default = "./examples/output/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.pdf",
    help = "output image file",
    metavar = "file"
  ),
  make_option(
    c("-r", "--region"),
    type = "character",
    default = "chr17:10958130-11017414",
    help = "region to plot in ucsc format",
    metavar = "region"
  ),
  make_option(
    c("-d", "--debug"),
    action = "store_true",
    default = FALSE,
    help = "debug mode, creates a tsv file of the alignments",
    dest = "debug"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
# create the output directory if it doesn't exist
dir.create(dirname(opt$output), showWarnings = FALSE)

processed = processRegion(opt$bam, opt$region)
savePlot(processed$g, opt$output)

if (opt$debug) {
  # write bamAll to a file named opt$output with a txt extension and the output extension removed
  g = gzfile(paste0(tools::file_path_sans_ext(opt$output), ".tsv.gz"),
             "w")
  write.table(
    processed$bamAll,
    file = g,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  close(g)
  
  # write adjustedDF to a file named opt$output with a txt extension and the output extension removed
  g = gzfile(paste0(tools::file_path_sans_ext(opt$output), ".adjusted.tsv.gz"),
             "w")
  write.table(
    processed$adjustedDF,
    file = g,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  close(g)
}
