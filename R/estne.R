#--------------------------------------
# Estimate effective populations size -
#--------------------------------------

checkSNP <- function(p0, pt, cov0, covt, truncAF=NA) {
  # return false if any of the following conditions is met: xi==0, xi==1
  # mask extreme allele frequencies if truncAF is unequal to 'NA'
  return(p0 != 0 & p0 != 1 & cov0 != 0 & covt != 0 & if(is.na(truncAF)) TRUE else p0 >= truncAF & p0 <= 1-truncAF)
}

estimateNe <- function(p0, pt, cov0, covt, t, ploidy=2, truncAF=NA, method="P.planI", Ncensus=NA, poolSize=rep(Ncensus, times=2), asList=FALSE) {

  # check if parameter 'method' has been set properly - stop execution otherwise
  mm <- match.arg(method, choices=c("P.planI", "P.planII",
                                    "JR.planI", "JR.planII",
                                    "W.planI", "W.planII",
                                    "P.alt.1step.planII", "P.alt.2step.planI", "P.alt.2step.planII"), several.ok=TRUE)
  if(length(mm) != length(method))
    stop("Unable to resolve the following method(s): ", paste("'", method[!method %in% mm], "'", sep="", collapse=", "))
  else
    method <- mm

  # check if poolSize parameter is set propperly - stop execution otherwise
  if(length(poolSize <- as.numeric(poolSize)) != 2)
    stop("Length of 'poolSize' parameter has to be equal to 2: length(poolSize) = ", length(poolSize))

  # remove SNPs for which Ne cannot/should not be estimated
  keep <- checkSNP(p0, pt, cov0, covt, truncAF=truncAF)
  xi <- p0[keep]
  yi <- pt[keep]
  zi = (xi + yi)/2
  n = length(xi)
  # coverage is divided by 2, because later 2*S0, 2*St correction term will be used
  s0 <- cov0[keep]
  st <- covt[keep]
  s0_mean <- mean(s0)
  st_mean <- mean(st)

  # estimate Ne using the specified method(s)
  res <- numeric(length=0)

  #--- Waples (1989) ---
  if(any(grepl("^W\\.plan(I|II)$", method))) {
    Fc <- ((xi-yi)^2)/(zi-xi*yi)
    # W_planI
    if(any(grepl("^W\\.planI$", method))) {
      Fc_planI <- (1/n)*sum( Fc - (1/s0 + 1/st) + 1/Ncensus )
      res <- c(res, Nw.planI=-t/(ploidy*log(1-Fc_planI)))
    }
    # W_planII
    if(any(grepl("^W\\.planII$", method))) {
      Fc_planII <- (1/n)*sum( Fc - (1/s0 + 1/st))
      res <- c(res, Nw.planII=-t/(ploidy*log(1-Fc_planII)))
    }
  }

  #--- Jorde and Ryman (2007) ---
  if(any(grepl("^JR\\.plan(I|II)$", method))) {
    F_nom <- (xi-yi)^2
    F_denom <- zi*(1-zi)
    n_harmonic <- 1/(1/s0 + 1/st)
    # JR_planI
    if(any(grepl("^JR\\.planI$", method))) {
      Fs_planI_nom <- sum(F_nom * (1 - 1/(4*n_harmonic) + 1/(4*Ncensus)) * n_harmonic * Ncensus + F_denom * (n_harmonic - Ncensus))/sum(F_denom * n_harmonic * Ncensus)
      Fs_planI_denom <- sum((4*F_denom + F_nom) * (1 - 1/st))/sum(4*F_denom)
      Fs_planI <- Fs_planI_nom/Fs_planI_denom
      res <- c(res, Njr.planI=-t/(ploidy*log(1-Fs_planI)))
    }
    # JR_planII
    if(any(grepl("^JR\\.planII$", method))) {
      Fs_planII_nom = sum(F_nom * (1-1/(4*n_harmonic)) * n_harmonic - F_denom )/sum( F_denom * n_harmonic )
      Fs_planII_denom = sum((4*F_denom + F_nom) * (1 - 1/st))/sum( 4*F_denom )
      Fs_planII = Fs_planII_nom/Fs_planII_denom
      res <- c(res, Njr.planII=-t/(ploidy*log(1-Fs_planII)))
    }
  }

  #--- Andreas Futschik ---
  if(any(grepl("^P\\.alt", method))) {
    # set correction term for 1-step sampling (plan II)
    if(any(grepl("^P\\.alt\\.1step\\.planII$", method))) {
      S <- sum(xi*(1-xi)*1/s0 + yi*(1-yi)*1/st)
      Ft <- (sum((xi-yi)^2) - S) / sum(xi*(1-xi))
      res <- c(res, Np.alt.1step.planII=-t/(ploidy*log(1-Ft)))
    }
    # ... 2-step sampling (plan II)
    if(any(grepl("^P\\.alt\\.2step\\.planII$", method))) {
      S <- sum(xi*(1-xi)*(1/s0+1/(ploidy*poolSize[1])-1/(s0*ploidy*poolSize[1])) + yi*(1-yi)*(1/st+1/(ploidy*poolSize[2])-1/(st*ploidy*poolSize[2])))
      Ft <- (sum((xi-yi)^2) - S) / sum(xi*(1-xi))
      res <- c(res, Np.alt.2step.planII=-t/(ploidy*log(1-Ft)))
    }
    # ... 2-step sampling (plan I)    ----------------TODO: NCENSUS AND POOLSIZE MUST BE MULTIPLIED WITH PLOIDY------------------------
    if(any(grepl("^P\\.alt\\.2step\\.planI$", method))) {
      S <- sum(xi*(1-xi)*(1/s0+1/(ploidy*poolSize[1])*(Ncensus-poolSize[1])/(Ncensus-1)-1/(s0*ploidy*poolSize[1])) + yi*(1-yi)*(1/st+1/(ploidy*poolSize[2])*(Ncensus-poolSize[2])/(Ncensus-1)-1/(st*ploidy*poolSize[2])))
      Ft <- (sum((xi-yi)^2) - S) / sum(xi*(1-xi))
      res <- c(res, Np.alt.2step.planI=-t/(ploidy*log(1-Ft)))
    }
  }

  #--- Agnes Jonas ---
  if(any(grepl("^P\\.plan", method))) {
    C0i <- 1/s0 + 1/(ploidy*poolSize[1]) - 1/(s0*ploidy*poolSize[1])
    Cti <- 1/st + 1/(ploidy*poolSize[2]) - 1/(st*ploidy*poolSize[2])

    # sampling plan II
    if(any(grepl("^P\\.planII$", method))) {
      Ft <- sum((xi - yi)^2 - (zi-xi*yi)*( C0i + Cti )) / sum( (zi-xi*yi) * (1 - Cti))
      res <- c(res, Np.planII=-t/(ploidy*log(1-Ft)))
    }
    # sampling plan I    ----------------TODO: NCENSUS MUST BE MULTIPLIED WITH PLOIDY------------------------
    if(any(grepl("^P\\.planI$", method))) {
      Ft <- sum((xi - yi)^2 * (1 - 1/(ploidy*Ncensus)) - (zi-xi*yi)*( C0i + Cti - 1/Ncensus)) / sum( (zi-xi*yi) * (1 - Cti))
      res <- c(res, Np.planI=-t/(ploidy*log(1-Ft)))
    }
  }

  # return either numeric vector or list
  if(asList)
    return(as.list(res))

  return(res)
}

estimateWndNe <- function(chr, pos, wndSize, p0, pt, cov0, covt, t, unit=c("bp", "SNP"), ploidy=2, truncAF=NA, method="P.planI", Ncensus=NA, poolSize=rep(Ncensus, times=2)) {

  # check unit parameter
  unit <- match.arg(unit)

  # organize all input vectors into one data table
  dataDt <- data.table(chr=as.character(chr), pos=pos, p0=p0, pt=pt, cov0=cov0, covt=covt)
  setkey(dataDt, chr)

  # remove SNPs for which Ne cannot/should not be estimated
  dataDt <- dataDt[checkSNP(p0, pt, cov0, covt),]

  # go through each chromosome and return results
  return(foreach(cc=unique(dataDt$chr), .combine=rbind) %do% {
    # get SNPs of current chromosome
    dataSubDt <- dataDt[cc]
    # assign window ID to each SNP based on specified windowSize and unit
    if(is.na(wndSize)) {
      dataSubDt[,wndID:=1]
    } else {
      dataSubDt[,wndID:=switch(unit, bp=floor(dataSubDt$pos/wndSize)+1, SNP=floor((rank(dataSubDt$pos)-1)/wndSize)+1)]
    }
    # estimate Ne for each window
    wndRes <- dataSubDt[,estimateNe(p0=p0, pt=pt, cov0=cov0, covt=covt, t=t, ploidy=ploidy, truncAF=truncAF, method=method, poolSize=poolSize, Ncensus=Ncensus, asList=TRUE), by=wndID]
    # add columh for chromosome name, window start/stop pos and number of SNPs per window
    wndRes <- wndRes[order(wndID)]
    wndMinMax <- switch(unit, bp=NA, SNP=dataSubDt[,list(min=min(pos), max=max(pos)),by=wndID])
    wndRes[,c("chr", "wndStart", "wndStop", "SNPs"):=list(cc,
                                                          switch(unit, bp=1+(wndID-1)*wndSize, SNP=wndMinMax$min[wndRes$wndID]),
                                                          switch(unit, bp=wndID*wndSize, SNP=wndMinMax$max[wndRes$wndID]),
                                                          as.integer(table(dataSubDt$wndID)[as.character(wndRes$wndID)]))]

    return(wndRes[,!grepl("^wndID$", colnames(wndRes)),with=FALSE])
  } )
}
