# for some reason crayon seems to install just fine in the docker image and can be loaded during image creation.
# But after the image is created, it can't be loaded. So, this installs it in case it's not already installed.
if (!require(crayon))
  install.packages("crayon")
library(crayon)
require(Rsamtools)
require(GenomicAlignments)


.unlist <- function (x)
{
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)) {
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

convertBam <- function(bam, param) {
  bam <- unname(bam) # names not useful in unlisted result
  elts <- setNames(bamWhat(param), bamWhat(param))
  lst <- lapply(elts, function(elt)
    .unlist(lapply(bam, "[[", elt)))
  lst = lst[lengths(lst) != 0]
  df <- do.call("data.frame", bam)
  return(df)
}

loadBam <- function(bamFile, param) {
  bam = scanBam(bamFile, param = param)
  bam = convertBam(bam = bam, param = param)
  return(bam)
}

parseAlignments <- function(bamFile, region) {
  which <-
    GRanges(seqnames  = region[1],
            ranges = IRanges(as.integer(region[2]), as.integer(region[3])))
  bamAll = loadBam(bamFile, param = ScanBamParam(what = scanBamWhat(), which = which))
  bamAll$flag = as.character(bamAll$flag)
  bamAll$sequenceLength = nchar(bamAll$seq)
  bamAll$cigarWidthAlongReferenceSpace = cigarWidthAlongReferenceSpace(bamAll$cigar)
  bamAll$cigarWidthAlongQuerySpace = cigarWidthAlongQuerySpace(bamAll$cigar)
  return(bamAll)
}

# bamAll$end = bamAll$pos + bamAll$swidth
# bamAll$proportionOfReadAlignedToRef = cigarWidthAlongQuerySpace(bamAll$cigar, after.soft.clipping = TRUE) / bamAll$sequenceLength
# bamAll$impliedEnd = bamAll$pos + bamAll$sequenceLength
# bamAll$diff = bamAll$impliedEnd - bamAll$end
# bamAll$op = explodeCigarOps(bamAll$cigar)
# bamAll$len = explodeCigarOpLengths(bamAll$cigar)
# bamAll$firstOp = sapply(bamAll$op, function(x)
#   x[1])
# bamAll$lastOp = sapply(bamAll$op, function(x)
#   x[length(x)])
# bamAll$firstLen = sapply(bamAll$len, function(x)
#   x[1])
# bamAll$lastLen = sapply(bamAll$len, function(x)
#   x[length(x)])
