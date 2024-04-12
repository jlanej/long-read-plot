require(GenomicAlignments)
require(Rsamtools)
require(ggplot2)
require(dplyr)
require(gridExtra)

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
  adjustedDF$line_type = "alignment-start-to-end"
  
  # create and incex for each unique qname
  sortedDF = adjustedDF[order(adjustedDF$pos),]
  indexDF = data.frame(qname = unique(sortedDF$qname),
                       qname_index = 1:length(unique(sortedDF$qname)))
  adjustedDF = inner_join(adjustedDF, indexDF, by = c("qname" = "qname"))
  adjustedDF$anonymousReadID = paste0("read_",
                                      adjustedDF$qname_index,
                                      "_alignment_",
                                      adjustedDF$alignment_number)
  
  
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
      rep(referenceSpaceInt, nrow(adjustedDF)),
      rep(readSpaceInt, nrow(adjustedDF)),
      rep(readSpaceInt, nrow(adjustedDF)),
      rep(referenceSpaceInt, nrow(adjustedDF))
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
getRearrangedDFStackRead <- function(adjustedDF) {
  rearranged = data.frame(
    pos = c(
      adjustedDF$adjustedPos,
      adjustedDF$adjustedPos,
      adjustedDF$adjustedPosEnd,
      adjustedDF$adjustedPosEnd
    ),
    type = c(
      adjustedDF$qname_index + rep(readSpaceInt + 1, nrow(adjustedDF)),
      adjustedDF$qname_index + rep(readSpaceInt, nrow(adjustedDF)),
      adjustedDF$qname_index + rep(readSpaceInt, nrow(adjustedDF)),
      adjustedDF$qname_index + rep(readSpaceInt + 1, nrow(adjustedDF))
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

getRearrangedDFStackRef <- function(adjustedDF) {
  rearranged = data.frame(
    pos = c(
      adjustedDF$pos,
      adjustedDF$pos,
      adjustedDF$end,
      adjustedDF$end
    ),
    type = c(
      rep(referenceSpaceInt - 1, nrow(adjustedDF)) - adjustedDF$qname_index,
      rep(referenceSpaceInt, nrow(adjustedDF)) - adjustedDF$qname_index,
      rep(referenceSpaceInt, nrow(adjustedDF)) - adjustedDF$qname_index,
      rep(referenceSpaceInt - 1, nrow(adjustedDF)) - adjustedDF$qname_index
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

readSpaceStakedInt = 105
readSpaceInt = 100
referenceSpaceInt = 50
referenceSpaceStackInt = 45


# TODO, add per read alignments on statcked portion, y-axis =readID
getParticlePlotStack <- function(gParticle,
                                 adjustedDF, alphaRibbons = 1) {
  stackedRead = getRearrangedDFStackRead(adjustedDF)
  stackedRef = getRearrangedDFStackRef(adjustedDF)
  g = gParticle + geom_polygon(
    data = stackedRead,
    aes(
      x = pos,
      y = type,
      group = qname,
      fill = alignment_number
    ),
    alpha = alphaRibbons
  )
  g = g + geom_polygon(
    data = stackedRef,
    aes(
      x = pos,
      y = type,
      group = qname,
      fill = alignment_number
      
    ),
    alpha = alphaRibbons
  )
  # g = g + scale_y_continuous(
  #   breaks = c(
  #     referenceSpaceStackInt,
  #     referenceSpaceInt,
  #     readSpaceInt,
  #     readSpaceStakedInt
  #   ),
  #   labels = c(
  #     "Reference Space Alignments",
  #     "Reference Space",
  #     "Read Space",
  #     "Read Space Alignments"
  #   )
  # )
  g = g + scale_y_continuous(
    breaks = c(referenceSpaceInt,
               readSpaceInt),
    labels = c("Reference Space",
               "Read Space")
  )
  
  return(g)
}

addCurves <- function(gParticle, adjustedDF,
                      curvature = 0.75) {
  # g = getParticlePlot(rearranged, rearrangedLines, adjustedDF, alphaRibbons, xlab)
  g = gParticle + geom_curve(
    data = adjustedDF,
    aes(
      x = adjustedPos ,
      y =  readSpaceInt,
      xend = adjustedPosEnd,
      yend = readSpaceInt,
      color = alignment_number,
      linetype = line_type
    ),
    alpha = 1,
    linewidth = 1,
    curvature = -1 * curvature,
    arrow = arrow()
  ) + geom_curve(
    data = adjustedDF,
    aes(
      x = pos ,
      y =  referenceSpaceInt,
      xend = end,
      yend = referenceSpaceInt,
      color = alignment_number,
      linetype = line_type
    ),
    alpha = 1,
    linewidth = 1,
    curvature = curvature,
    arrow = arrow()
  )
  g = g + scale_y_continuous(
    breaks = c(referenceSpaceInt,
               readSpaceInt),
    labels = c("Reference Space",
               "Read Space")
  ) + expand_limits(y = c(referenceSpaceStackInt - 15, readSpaceStakedInt +
                            15))
  return(g)
}
getParticlePlot <-
  function(rearranged,
           rearrangedLines,
           adjustedDF,
           alphaRibbons = 0.05,
           xlab = "Position") {
    g = ggplot(rearranged)
    g = g + geom_segment(
      data = rearrangedLines,
      aes(
        x = adjustedPos ,
        xend = refpos,
        y =  readSpaceInt,
        yend = referenceSpaceInt,
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
    
    g = g + scale_y_continuous(limits = c(referenceSpaceInt, readSpaceInt))
    g = g + xlab(xlab)
    # remove y axis label
    g = g + theme(axis.title.y = element_blank())
    g = g + scale_linetype_manual(values = c("twodash", "dotted", "solid"))
    return(g)
  }


getArrowPlot <-
  function(adjustedDF,
           linewidth,
           pointsize,
           yaxisFontSize) {
    adjustedDF$sameStart = abs(adjustedDF$adjustedPos - adjustedDF$pos) < 10
    g = ggplot(adjustedDF[which(!adjustedDF$sameStart),])
    geom_point(
      data = adjustedDF[which(adjustedDF$sameStart), ],
      aes(x = adjustedPos, y = anonymousReadID),
      color = "black",
      size = pointsize
    )
    g = g + geom_segment(
      data = adjustedDF,
      aes(
        x = pos,
        xend = end,
        y = anonymousReadID,
        yend = anonymousReadID,
        color = alignment_number
      ),
      alpha = 1,
      linewidth = 1,
      linetype = "dotted"
    )
    
    g = g + geom_point(
      data = adjustedDF,
      aes(x = end, y = anonymousReadID),
      color = "black",
      size = pointsize
    )
    #
    
    g = g + geom_segment(
      aes(
        x = pos,
        xend = adjustedPos,
        y = anonymousReadID,
        yend = anonymousReadID,
        color = alignment_number
      ),
      alpha = 1,
      linewidth = linewidth,
      arrow = arrow()
    )
    
    g = g + theme(panel.grid.minor = element_blank(),
                  panel.background = element_blank())
    
    g = g + theme(text = element_text(size = yaxisFontSize))
    
    return(g)
  }

# function to create an arrow plot for two sorting schemes
# 1 sort by the start of the alignment
# sort by the anonymous read ID
sortArrowPlots <- function(adjustedDF,
                           linewidth = .25,
                           pointsize = .5,
                           yaxisFontSize = 5) {
  base =  getArrowPlot(
    adjustedDF,
    linewidth = linewidth,
    pointsize = pointsize,
    yaxisFontSize = yaxisFontSize
  )
  gArrowStart = base + scale_y_discrete(limits = rev(adjustedDF[order(adjustedDF$qname_index), ]$anonymousReadID))
  gArrowReadID = base + scale_y_discrete(limits = rev(adjustedDF[order(adjustedDF$pos),]$anonymousReadID))
  results = list(gArrowStart = gArrowStart, gArrowReadID = gArrowReadID)
  return(results)
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
    gParticle = addCurves(getParticlePlot(rearranged, rearrangedLines, adjustedDF),
                          adjustedDF)
    gParticleStack = getParticlePlotStack(getParticlePlot(rearranged, rearrangedLines, adjustedDF),
                                          adjustedDF)
    gArrows = sortArrowPlots(adjustedDF)
    # stop()
    return(c(
      gArrows,
      list(
        gParticle = gParticle,
        gParticleStack = gParticleStack,
        bamAll = bamAll,
        adjustedDF = adjustedDF
      )
    ))
  }

savePlot <- function(l,
                     filename,
                     createIndividualPlots = TRUE,
                     width = 10,
                     singleheight =  5,
                     comboheight = 5,
                     dpi = "retina") {
  
  # ggsave("arrange2x2.pdf",,
  #        device = "pdf")
  
  if(createIndividualPlots) {
    invisible(mapply(
      ggsave,
      file = paste0(
        tools::file_path_sans_ext(filename),
        ".",
        names(l),
        ".",
        tools::file_ext(filename)
      ),
      plot = l,
      width = width,
      height = singleheight,
      dpi = dpi
    ))
  }
    
  ggsave(
    filename = filename,
    marrangeGrob(
      grobs = l,
      nrow = length(l),
      ncol = 1
    ),
    width = width,
    height = comboheight * length(l),
    dpi = dpi
  )
}