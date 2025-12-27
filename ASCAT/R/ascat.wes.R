runASCAT_WES <- function(PATH_BAF,
                         PATH_LOGR,
                         PATH_PON,
                         OUTDIR,
                         PENALTY=70,
                         PROP_HETERO = 0.05,
                         SAMPLENAME="Sample1",
                         GENDER="female")
{
    setwd(OUTDIR)
    require(ASCAT.sc)
    generateLogR.BAF_files(OUTDIR,
                           PATH_BAF,
                           PATH_LOGR,
                           PATH_PON)
    ## ##################################################
    ## STEP1 LOAD RAW DATA & CHANGE WORKDIR TO SAMPLEDIR
    ## ##################################################
    BAFin <- paste0(OUTDIR,"/BAFin.txt")
    LogRin <- paste0(OUTDIR,"/LogRin.txt")
    BAFM <- as.data.frame(data.table::fread(BAFin))
    rownames(BAFM) <- BAFM[,1]; BAFM <- BAFM[,-c(1)]
    LogRM <- as.data.frame(data.table::fread(LogRin))
    rownames(LogRM) <- LogRM[,1]; LogRM <- LogRM[,-c(1)]
    ## ##################################################
    ## ##################################################
    ## END STEP1 LOAD RAW DATA & CHANGE WORKDIR TO SAMPLEDIR
    ## ##################################################
    ## ##################################################
    ## STEP2 FIRST ASCAT RUN
    ## ##################################################
    ascat.bc  <-  ascat.loadData(paste0(OUTDIR,"/LogRin.txt"),paste0(OUTDIR,"/BAFin.txt"),
                                 NULL,
                                 NULL,
                                 gender=GENDER)
    gg <- ascat.predictGermlineGenotypes(ascat.bc, platform = "HTC_30X+")
    ## ##################################################
    ascat.bc <- ascat.aspcf(ascat.bc, ascat.gg=gg, penalty = PENALTY, out.prefix = "first")
    ascat.plot <- ascat.plotSegmentedData(ascat.bc, img.prefix = "first")
    ## ##################################################
    ascat.output <- ascat.runAscat(ascat.bc,
                                   img.prefix = "first",
                                   min_purity=.1,
                                   max_ploidy=2.9,
                                   gamma=1)
    if(is.null(ascat.output$segments))
        ascat.output <- ascat.runAscat(ascat.bc, img.prefix = "first",
                                       min_purity=.1, gamma=1,
                                       rho_manual = 1,## Change here
                                       psi_manual = 3.5) ## Change here
    ## ##################################################
    save(ascat.output, file="ASCAT.output.first.Rda")
    ## ##################################################
    save(ascat.bc, file="ASCAT.bc.first.Rda")
    ## ##################################################
    BAFPCfed <- as.data.frame(data.table::fread(dir(pattern=paste0(".BAF.PCFed.txt"))[1]))
    ## ##################################################
    ## END STEP2 FIRST ASCAT RUN
    ## ##################################################
    ## ##################################################
    ## STEP2b FIRST ASCAT.sc RUN
    ## ##################################################
    duplicates <- unlist(lapply(unique(LogRM[,1]),function(x) duplicated(LogRM[LogRM[,1]==x,2])))
    keep <- !is.na(LogRM[,3]) & !duplicates
    ## ##################################################
    track <- getTrackForAll.snp6.v2(LogRM[keep,3],
                                    CHRS=paste0("chr",LogRM[keep,1]),
                                    STARTS=LogRM[keep,2],
                                    ENDS=LogRM[keep,2]+1,
                                    allchr=c(1:22,"X"),
                                    ascat.output=ascat.output)
    ## ##################################################
    save(track,file="track.first.Rda")
    ## ##################################################
    GAMMA <- 1 ## 1 for sequencing - .55 for arrays
    ## ##################################################
    solution <- searchGrid(track,
                           purs=ascat.output$aberrantcellfraction+seq(-0.04, 0.1, 0.01),
                           ploidies=ascat.output$psi+seq(-0.3,0.3,0.01),
                           maxTumourPhi=7,
                           gamma=GAMMA)
    ## ##################################################
    save(solution,file="solution.first.Rda")
    ## ##################################################
    profile <- getProfile(fitProfile(track,
                                     purity=solution$purity,
                                     ploidy=solution$ploidy,
                                     gamma=GAMMA),
                          CHRS=c(1:22,"X"))
    ## ##################################################
    nprofile <- annotateProfiles(profile, BAFM, LogRM, gg)
    nprofile <- cbind(nprofile,
                      BAFconsider=considerBAF(nprofile[,"BAF"]))
    ## ##################################################
    ## END STEP2b FIRST ASCAT.sc RUN
    ## ##################################################
    ## ##################################################
    ## STEP3 RESCUE het SNPS in very imbalanced segments
    ## ##################################################
    newgg <- getmoreHetSNPs_where_depleted(gg, nprofile, BAFM, BAFPCfed, hetProp=PROP_HETERO)
    ## ##################################################
    ## END STEP3 RESCUE het SNPS in very imbalanced segments
    ## ##################################################
    ## ##################################################
    ## STEP4 SECOND ASCAT RUN -- play with min_purity and min/max_ploidy!
    ## ##################################################
    newgg$germlinegenotypes <- do.call("cbind",lapply(1,function(x) newgg$germlinegenotypes[,1]))
    colnames(newgg$germlinegenotypes) <- SAMPLENAME
    ascat.bc  <-  ascat.loadData(LogRin,BAFin, NULL, NULL, gender=rep(GENDER,1))
                                        #ascat.bc <- ascat.asmultipcf(ascat.bc, ascat.gg=newgg, penalty = PENALTY)
    ascat.bc <- ascat.aspcf(ascat.bc, ascat.gg=newgg, penalty = PENALTY)
    ascat.plot <- ascat.plotSegmentedData(ascat.bc, img.prefix = "second")
    ## ##################################################
    ascat.output <- ascat.runAscat(ascat.bc,
                                   img.prefix = "second",
                                   min_purity = .2,
                                   max_ploidy=2.7,
                                   gamma=GAMMA,
                                   write_segments = T)
    ## ##################################################
    ## ##################################################
    if(is.null(ascat.output$segments))
        ascat.output <- ascat.runAscat(ascat.bc, img.prefix = "second88",
                                       min_purity=.88,
                                       gamma=GAMMA,
                                       max_ploidy=8.5)
    if(is.null(ascat.output$segments))
        ascat.output <- ascat.runAscat(ascat.bc, img.prefix = "second14",
                                       min_purity=.99,
                                       gamma=GAMMA,
                                       rho_manual = 1,## Change here - sometimes slighlty <1 leads to best fit even for cell lines due to noise
                                       psi_manual = 4) ## Change here - play with this to find a working optima - visualise profiles
    ## ##################################################
    ## ##################################################
    save(ascat.output, file="ASCAT.output.second.Rda")
    ## ##################################################
    save(ascat.bc, file="ASCAT.bc.second.Rda")
    ## ##################################################
    ## ##################################################
    ## Write the output from the second ASCAT run to disk
    write.table(ascat.output$segments_raw,
                file=paste0(SAMPLENAME,"_ASCAT_final_output_raw.tsv"),
                sep="\t",col.names=T,row.names=F,quote=F)
    write.table(ascat.output$segments,
                file=paste0(SAMPLENAME,"_ASCAT_final_output.tsv"),
                sep="\t",col.names=T,row.names=F,quote=F)
    ## ##################################################
    ## END STEP4 SECOND ASCAT RUN
    ## ##################################################
}
