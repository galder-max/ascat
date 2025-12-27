## main function for rescuing SNPs in the very imbalanced segments
getmoreHetSNPs_where_depleted <- function(gg, nprofile, BAF, BAFPCfed, minprobes=100, hetProp=.2)
{
    require(GenomicRanges)
    rownames(BAFPCfed) <- BAFPCfed[,1]
    ##BAF <- BAF[grepl("SNP",rownames(BAF)),]
    grBAF <- GRanges(BAF[,1],IRanges(BAF[,2],BAF[,2]))
    grProf <- GRanges(nprofile[,1],
                      IRanges(as.numeric(as.character(nprofile[,2])),
                              as.numeric(as.character(nprofile[,3]))))
    ovs <- findOverlaps(grProf,grBAF)
    props <- nprofile$MeanProbesKept/nprofile$MeanProbes
    ##print(cbind(props,nprofile$MeanProbes))
    wwDepleted <- which(props<hetProp)
    kk <- try(median(props[nprofile[,"end"]-nprofile[,"start"]>20000000],na.rm=T),silent=T)
    if(inherits(kk,"try-error")) kk <- quantile(props,prob=.5,na.rm=T)
    hetProp <- max(hetProp,kk)
    okRetrieve <- !is.na(gg$germlinegenotypes[,1])
    for(ww in wwDepleted)
    {
        wwBAFinDepleted <- subjectHits(ovs)[queryHits(ovs)==ww]
        prop <- props[ww]
        medianBAF <- median(BAFPCfed[rownames(BAFPCfed)%in%rownames(BAF)[wwBAFinDepleted],2],na.rm=T)
        try({
            if(!is.null(medianBAF))
                if(!is.na(medianBAF))
                    if(nprofile$Nprobes[ww]>minprobes)
                        if(medianBAF>.75)
                        {
                            bbunsorted <- BAF[wwBAFinDepleted,3]
                            bbunsorted <- ifelse(bbunsorted<.5,1-bbunsorted,bbunsorted)
                            orderbb <- order(bbunsorted,decreasing=F)
                            extra <- wwBAFinDepleted[orderbb]
                            inds <- 1:round(length(wwBAFinDepleted)*(hetProp-prop))
                            keep <- rownames(BAF)[extra]%in%rownames(gg$germlinegenotypes)[okRetrieve]
                            rr <- rownames(BAF)[extra][keep][inds]
                            gg$germlinegenotypes[rownames(gg$germlinegenotypes)%in%rr,1] <- F
                        }
        })
    }
    gg
}
