getCorrectLogR <- function(lm, segs)
{
    require(GenomicRanges)
    gr1 <- GRanges(lm[,"Chromosome"],IRanges(lm[,"Position"],lm[,"Position"]))
    gr2 <- GRanges(seg[,"chr"],IRanges(seg[,"startpos"],seg[,"endpos"]))
    ovs <- findOverlaps(gr1,gr2)
    seg$meanlogr <- tapply(queryHits(ovs),subjectHits(ovs),function(x) median(lm[x,4]))[as.character(1:nrow(seg))]
    seg$sdlogr <- tapply(queryHits(ovs),subjectHits(ovs),function(x) mad(lm[x,4]))[as.character(1:nrow(seg))]
    mm <- lm[,4]-seg$meanlogr[subjectHits(ovs)]
    mm <- mm-mean(mm)
    mm/seg$sdlogr[subjectHits(ovs)]
}
