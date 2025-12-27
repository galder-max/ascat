generateLogR.BAF_files <- function(OUTDIR,
                                   PATH_BAF,
                                   PATH_LOGR,
                                   PATH_PON)
{
    load(PATH_PON)
    baf <- as.data.frame(data.table::fread(PATH_BAF))
    logr <- as.data.frame(data.table::fread(PATH_LOGR))
    allopos <- baf[,1]
    allopos <- allopos[allopos%in%rownames(logr.pon)]
    rownames(baf) <- baf[,1]
    rownames(logr) <- logr[,1]
    logr <- logr[allopos,]
    baf <- baf[allopos,]
    ## ##################################################
    treatmerge <- function(LOGRMERGE, allopos)
    {
        LOGRMERGE <- LOGRMERGE[allopos,]
        for(i in 1:ncol(LOGRMERGE))
        {
            med <- median(LOGRMERGE[,i],na.rm=T)
            LOGRMERGE[is.na(LOGRMERGE[,i]),i] <- med
        }
        LOGRMERGE
    }
    smoothNormals <- function(logr, logr.pon)
    {
        nisX = !grepl("^X_",rownames(logr.pon))
        notalreadyinpanel <- sapply(1:ncol(logr.pon),function(x) cor(logr.pon[nisX,x],logr[nisX])>.99)
        if(sum(notalreadyinpanel)>0) print(paste0(which(notalreadyinpanel)," already in panel"))
        cat(".")
        logr <- logr-median(logr)
        predicted <- lm(y ~ . - 1,
                        data = data.frame(y = logr,
                                          X = logr.pon[, !notalreadyinpanel])-1)$fitted.values
        logr <- logr-predicted
        logr <- logr-median(logr)
        logr
    }
    LOGRMERGE <- treatmerge(logr.pon,allopos)
    logr[,4] <- smoothNormals(logr[,4],LOGRMERGE)
    ##logr <- logr[,-c(1)]
    write.table(logr[,-c(1)],
                file=paste0(OUTDIR,"/LogRin.txt"),
                sep="\t",
                col.names=T,
                row.names=T,
                quote=F)
    write.table(baf[,-c(1)],
                file=paste0(OUTDIR,"/BAFin.txt"),
                sep="\t",
                col.names=T,
                row.names=T,
                quote=F)
}
## ##################################################
## END STEPINTERMED OVERSEG ASCAT RUN
## ##################################################

