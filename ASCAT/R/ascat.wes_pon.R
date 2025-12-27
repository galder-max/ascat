generate_PON_from_diploid_controls <- function(PATH_LOGR_CONTROLS_FEMALES,
                                               PATH_LOGR_CONTROLS_MALES=NULL,
                                               NMAX=50,
                                               FILEOUT = "pon.rda")
{
    paths_logr <- c(PATH_LOGR_CONTROLS_FEMALES,
                    PATH_LOGR_CONTROLS_MALES)
    logr <- lapply(paths_logr,function(x)
    {
        tt = as.data.frame(data.table::fread(x))
        if(x%in%PATH_LOGR_CONTROLS_MALES)
        {
            isX=grepl("^X_",tt[,1])
            tt[isX,4]=tt[isX,4]-mean(tt[isX,4],na.rm=T)
        }
        tt
    })
    allpos <- table(unlist(lapply(logr,function(x) x[,1])))
    allpos <- names(allpos[allpos==max(allpos)]) ##only keep positions present in all samples! check this is keeping most of them
    ##--------------------------------------------------------------------------
    ## order positions within ordered contigs
    allpos <- unlist(lapply(c(1:22,"X"),function(x) allpos[grepl(paste0("^",x,"_"),allpos)]))
    ##--------------------------------------------------------------------------
    logr.pon <- sapply(logr, function(x)
    {
        rownames(x)=x[,1]
        x=x[allpos,4]
        x-mean(x)
    })
    rownames(logr.pon) <- allpos
    ##--------------------------------------------------------------------------
    ## only run PCA if more than 50 samples
    ##--------------------------------------------------------------------------
    if(ncol(logr.pon)>NMAX){
        res<-prcomp(logr.pon)
        number.pcs <- if(ncol(logr.pon)>50) 50 else ncol(logr.pon)
        logr.pon<-res$x[,1:number.pcs]}
    ##--------------------------------------------------------------------------
    ## only run PCA if more than 50 samples
    ##--------------------------------------------------------------------------
    save(logr.pon, file=FILEOUT)
    return(logr.pon)
}
