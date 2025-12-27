## Function to use ASCAT segments for ASCAT.sc
getTrackForAll.snp6.v2 <- function (logr,
                                    CHRS,
                                    STARTS,
                                    ENDS,
                                    allchr = ALLCHR,
                                    ascat.output)
{
    require(GenomicRanges)
    segraw <- ascat.output$segments_raw
    ## smoothed <- my.smooth.CNA(logr, chrs=CHRS)
    smoothed <- logr
    lT <- lapply(allchr, function(chr) {
        cond <- CHRS == paste0("chr", chr)
        data.frame(start = STARTS[cond], end = ENDS[cond], records = logr[cond],
            fitted = logr[cond], smoothed = smoothed[cond])
    })
    lSe <- lapply(allchr, function(chr) {
        cond <- CHRS == paste0("chr", chr)
        list(starts = STARTS[cond], ends = ENDS[cond])
    })
    getSegs <- function(logr,chr,start,end, segraw)
    {
        if(chr==23) chr <- "X"
        keep1 <- gsub("chr","",segraw[,"chr"])==unique(gsub("chr","",chr))
        gr1 <- GRanges(gsub("chr","",segraw[keep1,"chr"]),
                       IRanges(segraw[keep1,"startpos"],
                               segraw[keep1,"endpos"]))
        gr2 <- GRanges(gsub("chr","",chr),
                       IRanges(start,end))
        ovs <- findOverlaps(gr1,gr2)
        data <- data.frame(chrom=chr,
                           maploc=start,
                           Sample.1=logr)
        output <- data.frame(ID="Sample.1",
                             chrom=chr,
                             loc.start=segraw[keep1,"startpos"],
                             loc.end=segraw[keep1,"endpos"],
                             num.mark=0,
                             seg.mean=NA)
        means <- tapply(1:length(ovs),queryHits(ovs),function(x) mean(logr[subjectHits(ovs)],na.rm=T))
        nmarks <- tapply(1:length(ovs),queryHits(ovs),function(x) length(x))
        output[as.numeric(names(means)),"num.mark"] <- nmarks
        output[as.numeric(names(means)),"seg.mean"] <- means
        list(data=data,
             output=output)
    }
    lSegs <- lapply(1:length(lT), function(x) {
        require(DNAcopy)
        cat(".")
        segments <- getSegs(lT[[x]]$smoothed, chr = paste0(x),
                             lSe[[x]]$starts, lSe[[x]]$ends, segraw)
    })
    names(lSegs) <- paste0(1:length(lT))
    tracks <- list(lCTS = lT, lSegs = lSegs)
    return(tracks)
}

