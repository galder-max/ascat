find_first_gaussian_peak <- function(x)
{
    xo = order(x,decreasing=F)[1:floor(length(x)/2)]
    xf=x[xo]
    keep=x<median(xf)*1.1+sd(xf)*1.6
    xf=x[keep]
    keep=x<median(xf)*1.1+sd(xf)*1.6
    first_peak_values=x[keep]
    return(first_peak_values)
}


smoothNormals <- function(logr, normals, L.more)
{
    notalreadyinpanel <- sapply(1:ncol(normals),function(x) cor(L.more[,x],logr)>.99)
    if(sum(notalreadyinpanel)>0) print(paste0(which(notalreadyinpanel)," already in panel"))
    cat(".")
    logr <- logr-median(logr)
    predicted <- lm(y ~ . - 1,
                    data = data.frame(y = logr,
                                      X = normals[, notalreadyinpanel])-1)$fitted.values
    logr <- logr-predicted
    logr <- logr-median(logr)
    logr
}

assignBinMeans <- function(values, reference_vector)
{
  # Ensure input vectors are of the same length
  if (length(values) != length(reference_vector)) {
    stop("Both input vectors must have the same length.")
  }
  # Compute quantile bins
  bins <- cut(values, breaks = quantile(values, probs = seq(0, 1, 0.1), na.rm = TRUE),
              include.lowest = TRUE, labels = FALSE)
  # Compute mean per bin
  bin_means <- tapply(reference_vector, bins, mean, na.rm = TRUE)
  # Assign bin means to the corresponding indices
  assigned_means <- bin_means[bins]
  return(assigned_means)
}


calculateMAPD <- function(x)
{
  # Ensure x is numeric
  if (!is.numeric(x)) {
    stop("Input must be a numeric vector.")
  }
  # Compute absolute pairwise differences
  pairwise_diffs <- abs(diff(x))
  # Compute median of absolute differences
  mapd_value <- median(pairwise_diffs, na.rm = TRUE)
  return(mapd_value)
}


calculateSDPD <- function(x) {
  # Ensure x is numeric
  if (!is.numeric(x)) {
    stop("Input must be a numeric vector.")
  }
  # Compute absolute pairwise differences
  pairwise_diffs <- abs(diff(x))
  # Compute standard deviation of absolute differences
  sd_value <- sd(pairwise_diffs, na.rm = TRUE)
  return(sd_value)
}

getCORRECTLOGR <- function(logr, logr.med)
{
    MEANLOGR <- runmed(logr,k=501)
    MEANLOGR.med <- apply(logr.med,1,median)
    CORRECTLOGR <- logr-(MEANLOGR.med-MEANLOGR)
    return(CORRECTLOGR)
    SDLOGR <- runmad(CORRECTLOGR,k=501)
    SDLOGR <- assignBinMeans(MEANLOGR, SDLOGR)
    CORRECTLOGR <- CORRECTLOGR/SDLOGR
}



## add BAF/LogR and density het SNPs to ASCAT profile
annotateProfiles <- function(profile, BAF, LogR, gg)
{
    library(GenomicRanges)
    grProf <- GRanges(profile[,1],IRanges(as.numeric(profile[,2]),
                                          as.numeric(profile[,3])))
    grBAF <-  GRanges(BAF[,1],IRanges(BAF[,2],BAF[,2]))
    ovs <- findOverlaps(grProf,grBAF)
    sizes <- as.numeric(profile[,"end"])-as.numeric(profile[,"start"])
    profile <- cbind(as.data.frame(profile),
                     sizes=sizes,
                     Nprobes=as.numeric(NA),
                     Nkept=as.numeric(NA),
                     MeanProbes=as.numeric(NA),
                     BAF=as.numeric(NA),
                     LogR=as.numeric(NA))
    BAFs <- tapply(1:length(ovs),queryHits(ovs),function(x)
    {
        keep <- 1:nrow(BAF)%in%subjectHits(ovs)[x]
        keep <- which(keep & gg$germlinegenotypes[,1])
        med <- median(BAF[keep,3],na.rm=T)
        return(med)
    })
    LogRs <- tapply(1:length(ovs),queryHits(ovs),function(x)
    {
        med <- median(LogR[subjectHits(ovs)[x],3],na.rm=T)
        med
    })
    Nprobes <- tapply(1:length(ovs),queryHits(ovs),function(x)
    {
        sum(!is.na(BAF[subjectHits(ovs)[x],3]))
    })
    Nkept <- tapply(1:length(ovs),queryHits(ovs),function(x)
    {
        sum(!is.na(BAF[subjectHits(ovs)[x],3]) & !gg$germlinegenotypes[subjectHits(ovs)[x],1],na.rm=T)
    })
    profile[as.numeric(names(BAFs)),"BAF"] <- BAFs
    profile[as.numeric(names(LogRs)),"LogR"] <- LogRs
    profile[as.numeric(names(Nprobes)),"Nprobes"] <- Nprobes
    profile[as.numeric(names(Nkept)),"Nkept"] <- Nkept
    profile[,"MeanProbes"] <- profile[,"Nprobes"]/sizes
    profile[,"MeanProbesKept"] <- profile[,"Nkept"]/sizes
    profile
}

## useful function for BAF prioritisation
considerBAF <- function(BAF, threshold=.2)
{
    sapply(BAF,function(x) if(is.na(x)) return(NA) else if(abs(x-.5)<threshold) return(x) else return(NA))
}


