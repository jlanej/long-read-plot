require(GenomicAlignments)
require(Rsamtools)
require(ggplot2)
require(dplyr)

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
  bamAll$cigarWidthAlongQuerySpaceAfterSoftClipping = cigarWidthAlongQuerySpace(bamAll$cigar, after.soft.clipping = TRUE)
  bamAll$end = bamAll$pos + bamAll$cigarWidthAlongReferenceSpace
  bamAll = addClipCounts(df = bamAll)
  bamAll <-
    bamAll %>% group_by(qname) %>% mutate(numAlignmentsForThisReadID = n())
  
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
  return(df)
}

# TODO make more R-like
getAdjustedDF <- function(df) {
  unique_qnames = unique(df$qname)
  adjustedDF = data.frame()
  df$adjustedPos = NA
  df$adjustedPosEnd = NA
  
  
  for (qname in unique_qnames) {
    dfSubset = df[df$qname == qname, ]
    if (nrow(dfSubset) > 1) {
      minSoftClip = min(dfSubset$LeftClipCount)
      minimumPosition = min(dfSubset[which(dfSubset$LeftClipCount == minSoftClip),]$pos)
      
      for (i in 1:nrow(dfSubset)) {
        sLens = explodeCigarOpLengths(dfSubset[i,]$cigar)
        sOps = explodeCigarOps(dfSubset[i,]$cigar)
        hasLeftSoft = FALSE
        if (sOps[[1]][[1]] == "S" || sOps[[1]][[1]] == "H") {
          dfSubset[i,]$adjustedPos = minimumPosition + sLens[[1]][[1]]
          dfSubset[i,]$adjustedPosEnd = dfSubset[i,]$adjustedPos + dfSubset[i,]$cigarWidthAlongReferenceSpace
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
  adjustedDF = adjustedDF[order(adjustedDF$LeftClipCount), ]
  adjustedDF$uniqueQname = make_unique(adjustedDF$qname, sep = "_aligment_#")
  adjustedDF$uniqueQname = ifelse(
    grepl("_aligment", adjustedDF$uniqueQname),
    adjustedDF$uniqueQname,
    paste0(adjustedDF$uniqueQname, "_aligment_#1")
  )
  adjustedDF$alignment_number = as.factor(gsub(".*_", "", adjustedDF$uniqueQname))
  adjustedDF$line_type="alignment-start-to-end"
  return(list(adjustedDF = adjustedDF, df = df))
}

make_unique <- function(vec, sep = "_") {
  vec[is.na(vec)] <- "NA"
  cs <- ave(vec == vec, vec, FUN = cumsum)
  vec[cs > 1] <- paste(vec[cs > 1], cs[cs > 1], sep = sep)
  vec
}

getRearrangedDF <- function(adjustedDF) {
  rearranged = data.frame(
    pos = c(
      adjustedDF$pos,
      adjustedDF$adjustedPos,
      adjustedDF$adjustedPosEnd,
      adjustedDF$end
    ),
    type = c(
      rep("reference-space", nrow(adjustedDF)),
      rep("read-space", nrow(adjustedDF)),
      rep("read-space", nrow(adjustedDF)),
      rep("reference-space", nrow(adjustedDF))
    ),
    qname = c(
      adjustedDF$uniqueQname,
      adjustedDF$uniqueQname,
      adjustedDF$uniqueQname,
      adjustedDF$uniqueQname
    )
  )
  rearranged$alignment_number = as.factor(gsub(".*_", "", rearranged$qname))
  return(rearranged)
}


getRearrangedLines <- function(adjustedDF) {
  rearrangedLines = data.frame(
    refpos = c(adjustedDF$pos,
               adjustedDF$end),
    adjustedPos = c(adjustedDF$adjustedPos,
                    adjustedDF$adjustedPosEnd),
    line_type = c(rep("alignment-start", nrow(adjustedDF)),
                  rep("alignment-end", nrow(adjustedDF))),
    qname = c(adjustedDF$uniqueQname,
              adjustedDF$uniqueQname)
  )
  rearrangedLines$alignment_number = as.factor(gsub(".*_", "", rearrangedLines$qname))
  return(rearrangedLines)
}

getParticlePlot <-
  function(rearranged,
           rearrangedLines,
           adjustedDF,
           alphaRibbons = 0.05,
           curvature = 0.75,
           xlab = "Position") {
    g = ggplot(rearranged)
    g = g + geom_segment(
      data = rearrangedLines,
      aes(
        x = adjustedPos ,
        xend = refpos,
        y =  "read-space",
        yend = "reference-space",
        color = alignment_number,
        linetype = line_type
      ),
      alpha = 1,
      linewidth = 1
    ) + geom_polygon(aes(
      x = pos,
      y = type,
      group = qname,
      fill = alignment_number
    ),
    alpha = alphaRibbons)
    g = g + guides(colour = guide_legend(override.aes = list(alpha = 1)))
    
    g = g + theme(panel.grid.minor = element_blank(),
                  panel.background = element_blank())
    
    g = g + geom_curve(
      data = adjustedDF,
      aes(
        x = adjustedPos ,
        y =  "read-space",
        xend = adjustedPosEnd,
        yend = "read-space",
        color = alignment_number,linetype = line_type
      ),
      alpha = 1,
      linewidth = 1,
      curvature = -1 * curvature,
      arrow = arrow()
    ) + geom_curve(
      data = adjustedDF,
      aes(
        x = pos ,
        y =  "reference-space",
        xend = end,
        yend = "reference-space",
        color = alignment_number,linetype = line_type
      ),
      alpha = 1,
      linewidth = 1,
      curvature = curvature,
      arrow = arrow()
    )
    g = g + scale_y_discrete(limits = c("reference-space", "read-space"))
    g = g + xlab(xlab)
    # remove y axis label
    g = g + theme(axis.title.y = element_blank())
    g=g+scale_linetype_manual(values=c("twodash","dotted","solid"))
    return(g)
  }

processRegion <-
  function(bamFile,
           ucscRegion = "chr1:1-1000",
           minAlignments = 2) {
    region = strsplit(ucscRegion, ":|-")[[1]]
    bamAll = parseAlignments(bamFile, region)
    bamAll = bamAll[which(bamAll$numAlignmentsForThisReadID >= minAlignments), ]
    adjusted = getAdjustedDF(df = bamAll)
    adjustedDF = adjusted$adjustedDF
    bamAll = adjusted$df
    rearranged = getRearrangedDF(adjustedDF)
    rearrangedLines = getRearrangedLines(adjustedDF)
    gParticle = getParticlePlot(rearranged, rearrangedLines, adjustedDF)
    return(list(
      gParticle = gParticle,
      bamAll = bamAll,
      adjustedDF = adjustedDF
    ))
  }

savePlot <- function(g,
                     filename,
                     width = 10,
                     height = 5,
                     dpi = "retina") {
  ggsave(
    filename = filename,
    plot = g,
    width = width,
    height = height,
    dpi = dpi
  )
}