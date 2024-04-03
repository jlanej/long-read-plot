require(GenomicAlignments)
require(Rsamtools)
require(ggplot2)

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
  bamAll$cigarWidthAlongReferenceSpaceAfterSoftClipping = cigarWidthAlongReferenceSpace(bamAll$cigar)
  bamAll$cigarWidthAlongQuerySpaceAfterSoftClipping = cigarWidthAlongQuerySpace(bamAll$cigar, after.soft.clipping = TRUE)
  addClipCounts(df = bamAll)
  return(bamAll)
}

addClipCounts <- function(df) {
  df$LeftClipCount = 0
  df$RightClipCount = 0
  for (i in 1:nrow(df)) {
    sLens = explodeCigarOpLengths(df[i, ]$cigar)
    sOps = explodeCigarOps(df[i, ]$cigar)
    if (sOps[[1]][[1]] == "S" | sOps[[1]][[1]] == "H") {
      df[i, ]$LeftClipCount = sLens[[1]][[1]]
    }
    if (sOps[[1]][[length(sOps[[1]])]] == "S" ||
        sOps[[1]][[length(sOps[[1]])]] == "H") {
      df[i, ]$RightClipCount = sLens[[1]][[length(sOps[[1]])]]
    }
  }
}

getAdjustedDF <- function(df) {
  unique_qnames = unique(df$qname)
  adjustedDF = data.frame()
  for (qname in unique_qnames) {
    dfSubset = df[df$qname == qname, ]
    if (nrow(dfSubset) > 1) {
      minSoftClip = min(dfSubset$LeftSoftClipCount)
      minimumPosition = min(dfSubset[which(dfSubset$LeftSoftClipCount == minSoftClip),]$pos)
      
      for (i in 1:nrow(dfSubset)) {
        sLens = explodeCigarOpLengths(dfSubset[i,]$cigar)
        sOps = explodeCigarOps(dfSubset[i,]$cigar)
        hasLeftSoft = FALSE
        if (sOps[[1]][[1]] == "S" || sOps[[1]][[1]] == "H") {
          dfSubset[i,]$adjustedPos = minimumPosition + sLens[[1]][[1]]
          dfSubset[i,]$adjustedPosEnd = dfSubset[i,]$adjustedPos + dfSubset[i,]$swidth
          hasLeftSoft = TRUE
        }
        if (!hasLeftSoft) {
          dfSubset[i,]$adjustedPos = dfSubset[i,]$pos
          dfSubset[i,]$adjustedPosEnd = dfSubset[i,]$end
        }
        adjustedDF = rbind(adjustedDF, dfSubset[i,])
      }
    }
  }
  return(adjustedDF)
}