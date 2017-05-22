# ----------------------------------
# Methods for estimate s/h results -
# ----------------------------------

print.estsh <- function(x, ...) {
  if(!is(x, class="estsh"))
    stop("'x' is not an object of class 'estsh'")

  hfix <- !is.null(x$h.given)

  # title
  if(x$used == "LLS") {
    cat("\n\tEstimation of s and p0 with linear least squares\n\n")
  } else if(x$used == "NLS") {
    if(hfix) {
      cat("\n\tEstimation of s with nonlinear least squares\n\n")
    } else {
      cat("\n\tEstimation of s and h with nonlinear least squares\n\n")
    }
  }

  # parameters
  cat("Ne:", x$Ne, if(x$haploid) "haploid" else "diploid", "individuals\n", sep=" ")
  cat(" t: ", paste(x$t, sep="", collapse="-"), "\n", sep="")
  if(!is.null(x$cov))
    cat("cov: ", paste(round(quantile(x$cov, probs=c(0.25, 0.75))), sep="", collapse="-"), " (IQR)\n", sep="")
  if(hfix)
    cat(" h: ", x$h.given, "\n",sep="")
  cat("\n")


  # results
  cat("s = ", x$s, sep="")
  if(!is.null(x$p0))
    cat(", p0 = ", x$p0, sep="")
  if(x$used == "NLS" && !is.null(x$h))
    cat(", h = ", x$h, sep="")
  if(!is.null(x$p.value))
    cat(", p.value = ", if(!is.na(x$p.value) && x$p.value*x$N.pval < 10) paste0("< ", 1/x$N.pval*10) else x$p.value, sep="")

  cat("\n")
  invisible(x)
}

coef.estsh <- function(object, ...) {
  if(!is(object, class="estsh"))
    stop("'object' is not an object of class 'estsh'")

  # combine all estimated parameters
  out <- c(s=object$s)
  if(!is.null(h <- object$h))
    out <- c(out, h=h)
  if(!is.null(p0 <- object$p0))
    out <- c(out, p0=p0)

  return(out)
}

confint.estsh <- function(object, parm, level = 0.95, N.ci = 1000, warn = FALSE, ...) {
  if(!is(object, class="estsh"))
    stop("'object' is not an object of class 'estsh'")

  # get coefficients (names and values)
  cf <- coef(object)
  pnames <- names(cf)
  # determine for which ones CIs should be returned
  if (missing(parm)) {
    parm <- pnames
  } else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }

  # simulate AF trajectories with parameters specified in 'object'
  repl <- nrow(object$traj)
  p0 <- if(!is.null(object$p0)) object$p0 else object$ctraj$af[1]
  Ne <- object$Ne
  t <- object$t
  haploid <- object$haploid
  s <- cf["s"]
  hfix <- !is.null(object$h.given)
  h <- if(hfix) object$h.given else object$h
  simTraj <- wf.traj(p0=rep(p0, times=repl*N.ci), Ne=Ne, t=t, s=s, h=h, haploid=haploid)
  # add sampling noise according to 'object$cov'
  if(!is.null(object$cov)) {
    smpl <- sample.alleles(p=as.vector(simTraj), size=as.vector(object$cov), mode="coverage", ploidy=if(haploid) 1 else 2)
    simTraj <- matrix(smpl, ncol=ncol(simTraj), dimnames=dimnames(simTraj))
  }

  # compute bias of consensus trajectory
  cBias <- fixationBias(traj=object$traj, Ne=Ne, t=t, N.sim=object$N.ctraj)

  method <- object$method
  sim.param <- list(s=numeric(N.ci))
  if(method == "LLS")
    sim.param$p0 = numeric(N.ci)
  if(!hfix)
    sim.param$h = numeric(N.ci)

  # estimate s/h for each simulated locus
  for(l in 1:N.ci) {
    # compute consensus trajectory
    cT <- consensus.traj(traj=simTraj[seq((l-1)*repl+1, l*repl),], t=t, cov=if(!is.null(object$cov)) object$cov else NA, bias=cBias)
    # determine if LLS or NLS should be used
    if(method == "LLS" || (method == "automatic" && useLLS(ctraj=cT, h=h, haploid=haploid))) {
      # perofrm linear regression (approximation)
      param <- lm.s(ctraj=cT, haploid=haploid)
      sim.param$s[l] <- param["s"]
      sim.param$p0[l] <- param["p0"]
    } else {
      # perform non-linear least-squares regression
      param <- nls.sh(ctraj=cT, Ne=Ne, h=if(hfix) h else NA, haploid=haploid, s.start=0.1, h.start=0.5)
      sim.param$s[l] <- param["s"]
      if(!hfix)
        sim.param$h[l] <- param["h"]
    }
  }

  # compute quantile thresholds
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- paste(format(100 * a, trim=TRUE, scientific=FALSE, digits=3), "%")
  ci <- matrix(NA, nrow=length(parm), ncol=2, dimnames=list(parm, pct))

  if(sum(is.na(sim.param$s)) > max(N.ci*(1-level), 0.05)) {
    if(warn)
      warning("Confidence interval might be inaccurate under the given scenario.")
    else
      return(ci)
  }

  # get sample from distribution of d = x* - x
  for(p in parm) {
    x <- cf[p]
    dq <- quantile(sim.param[[p]] - x, a, na.rm=TRUE)
    ci[p,] <- c(x - dq[2], x - dq[1])
  }

  return(ci)
}

print.compsh <- function(x, ...) {
  if(!is(x, class="compsh"))
    stop("'x' is not an object of class 'compsh'")

  # title
  hfix <- !is.null(x$h.given)
  if(hfix) {
    cat("\n\tCompare estimates of s between replicates\n\n")
  } else {
    cat("\n\tCompare estimates of s and h between replicates\n\n")
  }

  # parameters
  cat(" Ne:", x$Ne, if(x$haploid) "haploid" else "diploid", "individuals\n", sep=" ")
  cat("  t: ", paste(x$t, sep="", collapse="-"), "\n", sep="")
  if(!is.null(x$cov))
    cat("cov: ", paste(quantile(x$cov, probs=c(0.25, 0.75)), sep="", collapse="-"), " (IQR)\n", sep="")
  if(hfix)
    cat(" h: ", x$h.given, "\n",sep="")
  cat("\n")


  # results
  cat("s = ", x$s, ", var(s) = ", var(x$s.sep),
      ", p.value = ", if(x$s.p.value * x$N.pval < 10) paste0("< ", 1/x$N.pval*10) else x$s.p.value, "\n", sep="")
  if(!is.null(x$h)) {
    cat("h = ", x$h, ", var(h) = ", var(x$h.sep),
        ", p.value = ", if(x$h.p.value * x$N.pval < 10) paste0("< ", 1/x$N.pval*10) else x$h.p.value, "\n", sep="")
  }

  invisible(x)
}

# -----------------------------------------------------
# Estimate s and h for temporal allele frequency data -
# -----------------------------------------------------

scaleAF <- function(af, method=c("logit", "asin")) {
  # scale allele frequencies according to method
  method <- match.arg(method)
  switch(method,
         logit = log(af / (1-af)),
         asin = asin(sqrt(af)))
}

fixationBias <- function(traj, t, Ne, N.sim=10000) {

  if(length(traj) == 0)
    stop("No trajectory provided.")

  # if only 1 replicate is provided then return '0'
  if(is.vector(traj)) {
    return(0)
  } else {
    traj <- as.matrix(traj)
    if(nrow(traj) == 1)
      return(0)
  }

  # if no simulations should be performed then return '0'
  if(is.na(N.sim) || N.sim < 1)
    return(0)

  if(ncol(traj) != length(t <- as.numeric(t)))
    stop("Number of columns in 'traj' has to be equal to the length of 't'.")

  if(length(Ne <- as.numeric(Ne)) != 1 || length(N.sim <- as.numeric(N.sim)) != 1)
    stop("Both 'Ne' and 'N.sim' have to be of length 1.")

  p.mean <- colMeans(traj)

  if(p.mean[1] == 0)
    return(rep.int(0, length(p.mean)))

  # compute mean AF conditioning on fixation
  traj[rowCummaxs(traj[,ncol(traj):1])[,ncol(traj):1] == 0] <- NA
  p.mean.fix <- ifelse(colAlls(is.na(traj)), 0, colMeans(traj, na.rm=TRUE))

  t.diff <- diff(t)
  pt.mean <- p.mean[1]
  pt <- rep(p.mean[1], times=N.sim <- as.numeric(N.sim))
  bias <- rep.int(0, length(p.mean))

  # for all except last mean AF
  for(i in seq(1, length(p.mean)-1)) {
    # if AF equals '0' then no bias correction is required
    if(length(pt) == 0) {
      bias[i+1] <- bias[i]
      pt <- rep(p.mean[i+1], times=N.sim)
    # if AF equals '1' then no bias correction is required
    #} else if(p.mean.fix[i] == 1) {
    #  bias[i+1] <- 0
    # otherwise assess bias of conditioning on fixation
    } else {
      # simulate AF after time interval of random drift
      pt <- wf.traj(p0=pt, Ne=Ne, t=t.diff[i], s=0.0)
      # compute bias of consensus trajectory
      pt <- pt[pt != 0]
      bias[i+1] <-  if(length(pt) == 0) bias[i] else bias[i] + mean(pt) - pt.mean

      # update AF according to observed mean AF
      pt <- pt + p.mean.fix[i+1] - mean(pt)#p.mean.fix[i] - bias[i+1] + bias[i]
      pt <- ifelse(pt < 0, 0, ifelse(pt > 1, 1, pt))
      pt <- pt[pt > 0]
      pt.mean <- mean(pt)
    }
  }

  return(bias)
}

consensus.traj <- function(traj, t, cov, Ne, N.sim=10000, bias=NA) {

  if(length(traj) == 0)
    stop("No trajectory provided.")

  if(is.vector(traj)) {
    traj <- matrix(traj, nrow=1)
  } else {
    traj <- as.matrix(traj)
  }

  if(is.vector(cov)) {
    cov <- matrix(cov, nrow=1)
  } else {
    cov<- as.matrix(cov)
  }

  if(ncol(traj) != length(t <- as.numeric(t)))
    stop("Number of columns in 'traj' has to be equal to the length of 't'.")

  if(!missing(cov) && !all(is.na(cov)) && ncol(cov) != length(t))
    stop("Number of columns in 'cov' has to be equal to the length of 't'.")

  if((!missing(Ne) && length(Ne <- as.numeric(Ne)) != 1) || length(N.sim <- as.numeric(N.sim)) != 1)
    stop("Both 'Ne' and 'N.sim' have to be of length 1.")

  ctraj.out <- NULL
  # if there is only a single replicate then write it to output object
  if(nrow(traj) == 1) {
    ctraj.out <- list(af=traj[,order(t)], t=t)
    # if specified than provide coverage values in result object
    if(!missing(cov) && !all(is.na(cov))) {
      ctraj.out$cov <- cov[,order(t)]
    }

  # if there are multiple replicates then combine them
  } else {
    # order AFs increasing by generations
    traj <- traj[,order(t)]
    p0 <- mean(traj[,1])

    # get mean trajectory ignoring replicates, where the allele is lost for all successive time points
    # -> if allele is lost in all replicates of a time point then '0' is returned
    mask <- rowCummaxs(traj[,ncol(traj):1])[,ncol(traj):1] == 0
    traj[mask] <- NA
    ct <- ifelse(colAlls(is.na(traj)), 0, colMeans(traj, na.rm=TRUE))

    cc <- NULL
    # if coverage values are specified then apply masking procedure
    if(!missing(cov) && !all(is.na(cov))) {
      cov <- cov[,order(t)]
      cov[mask] <- NA
      cc <- ifelse(colAlls(is.na(cov)), 0, colSums(cov, na.rm=TRUE))
    }

    # remove bias that results from conditioning on fixation and ensure that 0 <= af <= 1
    ct <- ct - if(all(is.na(bias))) fixationBias(traj=traj, t=t, Ne=Ne, N.sim=N.sim) else bias
    ct <- ifelse(ct < 0, 0, ifelse(ct > 1, 1, ct))
    ct[1] <- p0

    # add ctraj and coverage to output object
    ctraj.out <- list(af=ct, t=t)
    if(!is.null(cc))
      ctraj.out$cov <- cc
  }

  class(ctraj.out) <- "ctraj"
  return(ctraj.out)
}

useLLS <- function(ctraj, h, haploid, p.min=0.10) {

  if(!haploid && !is.na(h) && h != 0.5)
    return(FALSE)

  if(!inherits(ctraj, "ctraj") || length(ctraj$af) == 0 || length(ctraj$t) == 0)
    stop("The parameter 'ctraj' is not specified properly. It has to be generated with consensus.traj and provide non-zero length.")

  if(length(p.min <- as.numeric(p.min)) != 1)
    stop("Both 'p.min' has to be of length 1.")

  # allele frequencies of '1' and '0' will produce +Inf/-Inf so we need to get rid of those
  ikeep <- which(ctraj$af != 0 & ctraj$af != 1)

  # if locus is polymorphic at less than 3 time points
  if(length(ikeep) < 3) {
    return(TRUE)
  }

  # scale allele frequencies and fit linear model
  fit <- lm(scaleAF(ctraj$af[ikeep], method="logit") ~ poly(ctraj$t[ikeep], degree=2, raw=TRUE))
  # compute p-value of quadratic term
  p <- summary(fit)$coefficients[3,4]

  # if p-value is smaller than 'p.min'
  if(!is.na(p) && p < p.min)
    return(FALSE)
  else
    return(TRUE)
}

nls.sh <- function(ctraj, Ne, h, haploid, s.start=0.1, h.start=0.5, approximate=FALSE) {
  af <- ctraj$af
  # obtaine initial NLS fit ignoring any warning messages
  fit <- tryCatch(suppressWarnings( {
    if(is.na(h)) {
      nls(af ~ wf.traj(p0=af[1], Ne=NA, t=ctraj$t, s=sEst, h=hEst, haploid=haploid, approximate=approximate), start=list(sEst=s.start, hEst=h.start), control=nls.control(warnOnly=TRUE))
    } else {
      nls(af ~ wf.traj(p0=af[1], Ne=NA, t=ctraj$t, s=sEst, h=h, haploid=haploid, approximate=approximate), start=list(sEst=s.start), control=nls.control(warnOnly=TRUE))
    }
  } ),
  error=function(e) {
    NULL
  } )

  res <- if(is.na(h)) c(s=NA_real_, h=NA_real_) else c(s=NA_real_)

  s.lim <- c(-Inf, +Inf)
  h.lim <- c(-Inf, +Inf)

  # if no error occurred in NLS
  if(!is.null(fit)) {
    # get estimates
    param <- coef(fit)

    s.fix <- NA
    h.fix <- NA
    # if s-estimate is outside of limits
    if(param["sEst"] < s.lim[1] || param["sEst"] > s.lim[2]) {
      s.fix <- if(param["sEst"] < s.lim[1]) s.lim[1] else s.lim[2]
    }
    # if h-estimate is outside of limits
    if(is.na(h) && (param["hEst"] < h.lim[1] || param["hEst"] > h.lim[2])) {
      h.fix <- if(param["hEst"] < h.lim[1]) h.lim[1] else h.lim[2]
    }

    # if both s/h should be estimated but only one is outside of limits
    if(is.na(h) && is.na(s.fix) + is.na(h.fix) == 1) {
      # re-estimate s/h using one fixed parameter value
      fit <- tryCatch(suppressWarnings( {
        fit <- if(!is.na(s.fix)) {
          nls(ctraj ~ wf.traj(p0=ctraj[1], Ne=NA, t=t, s=s.fix, h=hEst, haploid=haploid), start=list(hEst=h.start), control=nls.control(warnOnly=TRUE))
        } else {
          nls(ctraj ~ wf.traj(p0=ctraj[1], Ne=NA, t=t, s=sEst, h=h.fix, haploid=haploid), start=list(sEst=s.start), control=nls.control(warnOnly=TRUE))
        }
      } ),
      error=function(e) {
        NULL
      } )

      # if no error occurred in NLS
      if(!is.null(fit)) {
        # overwrite old estimates
        param <- coef(fit)
      }
    }

    # provide estimates only if NLS converged successfully
    if(!is.null(fit) && fit$convInfo$stopCode == 0) {
      res["s"] <- if(is.na(s.fix)) param["sEst"] else s.fix
      res["s"] <- if(res["s"] < s.lim[1]) s.lim[1] else if(res["s"] > s.lim[2]) s.lim[2] else res["s"]
      if(is.na(h)) {
        res["h"] <- if(is.na(h.fix)) param["hEst"] else h.fix
        res["h"] <- if(res["h"] < h.lim[1]) h.lim[1] else if(res["h"] > h.lim[2]) h.lim[2] else res["h"]
      }
    }
  }

  return(res)
}

lm.s <- function(ctraj, haploid, maxiter=10, tol=0.01) {

  # make sure that 'ctraj' is set properly
  if(!inherits(ctraj, "ctraj") || length(ctraj$af) == 0 || length(ctraj$t) == 0)
    stop("The parameter 'ctraj' is not specified properly. It has to be generated with consensus.traj and provide non-zero length.")

  # mask uninformative time points
  caf <- ctraj$af
  ct <- ctraj$t
  caf.diff <- abs(diff(caf))
  caf.maxdiff <- max(caf.diff, na.rm=TRUE)
  caf.rle <- sapply(rle(caf.diff < caf.maxdiff*0.01), tail, n=1, USE.NAMES=FALSE)
  if(!is.na(caf.rle["values"]) && caf.rle["values"] == 1) {
    irem <- seq(length(caf)-caf.rle["lengths"]+1, length(caf))
    caf <- caf[-irem]
    ct <- ct[-irem]
  }

  # allele frequencies of '1' and '0' will produce +Inf/-Inf so we need to get rid of those
  ikeep <- which(caf != 0 & caf != 1)

  # if all time points were masked
  if(length(ikeep) < 2) {
    return(c(s=NA_real_, p0=NA_real_))
  }
  caf <- caf[ikeep]
  ct <- ct[ikeep]

  # scaled allele frequencies, fit linear model and estimate s/p0
  lm.res <- lm(scaleAF(caf, method="logit") ~ ct)
  s <- unname(if(haploid) coef(lm.res)[2] else coef(lm.res)[2]*2)
  p0 <- unname(1/(exp(-coef(lm.res)[1]) + 1))

  # if possible then correct bias for large s owing to continuous time approximation
  if(s != 0 && p0 != 0 && p0 != 1 && !is.na(maxiter) && maxiter >= 1) {

    i <- 0
    s.delta <- 0
    s.est <- NA_real_
    # iteratively approximate bias
    repeat {
      # simulate trajectory for given estimate
      trj.sim <- wf.traj(p0=p0, Ne=NA, t=ct, s=s+s.delta, haploid=haploid)
      # stop if locus is polymorphic at less than 2 time points
      if(sum(trj.sim > 0 & trj.sim < 1) < 2) {
        warning("Correction not successfull. Selection estimates may be inaccurate.")
        break
      }
      # re-estimate parameters from simulations
      fit <- lm(scaleAF(trj.sim, method="logit") ~ ct)
      # compute new parameters from deviation between simulated and estimated s
      s.est <- unname(if(haploid) coef(fit)[2] else coef(fit)[2]*2)
      s.delta <- s + s.delta - s.est
      i <- i + 1

      # stop if either max. number of iterations or tolerance criterion is reached
      if(i >= maxiter || is.na(s.est) || abs((s.est-s)/s) <= tol)
        break
    }
    if(i >= maxiter && abs((s.est-s)/s) > tol)
      warning("Reached max. number of iterations without fullfilling convergence criterion. Selection estimates may be inaccurate.")

    # correct s-estimate
    s <- s + s.delta
  }

  return(c(s=s, p0=p0))
}

estimateSH <- function(traj, t, Ne, haploid=FALSE, h=NA, N.ctraj=0, simulate.p.value=FALSE, N.pval=1000, cov=NA, approximate=TRUE, method=c("LLS", "NLS", "automatic")) {
  if(is.vector(traj)) {
    traj <- matrix(traj, nrow=1)
  } else {
    traj <- as.matrix(traj)
  }
  if(ncol(traj) <= 1)
    stop("Allele frequencies at more than one time point have to be specified.")
  if(ncol(traj) != length(t <- as.numeric(t)))
    stop("Number of columns in 'traj' has to be equal to the length of 't'.")
  if(length(Ne <- as.numeric(Ne)) != 1)
    stop("Length of 'Ne' has to be equal to 1.")

  method <- match.arg(method)
  if((is.na(h) || h != 0.5) && method == "LLS")
    stop("If LLS is used to estimate s, then h has to equal 0.5.")

  # compute the consensus trajectory
  cBias <- fixationBias(traj=traj, t=t, Ne=Ne, N.sim=N.ctraj)
  cTraj <- consensus.traj(traj=traj, t=t, cov=cov, bias=cBias)
  p0 <- cTraj$af[1]

  param <- NULL
  used <- ""
  if(method == "LLS" || (method == "automatic" && useLLS(ctraj=cTraj, h=h, haploid=haploid))) {
    # perofrm linear regression (approximation)
    param <- lm.s(ctraj=cTraj, haploid=haploid, maxiter=if(approximate) 0 else 10)
    used <- "LLS"
  } else {
    # perform non-linear least-squares regression
    param <- nls.sh(ctraj=cTraj, Ne=Ne, h=h, haploid=haploid, s.start=0.1, h.start=0.5, approximate=approximate)
    used <- "NLS"
  }

  # generate result
  estsh.out <- list(traj=traj, t=t, Ne=Ne, haploid=haploid, ctraj=cTraj, N.ctraj=N.ctraj, cov=if(all(is.na(cov))) NULL else cov, method=method, used=used, s=unname(param["s"]))
  if("p0" %in% names(param)) {
    estsh.out$p0 <- unname(param["p0"])
  }
  if(is.na(h)) {
    if(used == "LLS")
      estsh.out$h.given <- 0.5
    else
      estsh.out$h <- unname(param["h"])
  } else {
    estsh.out$h.given <- h
  }

  # estimate p-value under the nullhypothesis of s=0
  if(simulate.p.value) {
    # make sure h is specified
    if(is.na(h)) {
      stop("A p-value can only be estimated if 'h' is specified and not 'NA'.")
    }

    # set p-value to 'NA' by default
    pval <- NA
    # if not'NA' then provide coverage in output
    if(any(!is.na(cov))) {
      # make sure that 'cov' is set properly
      if(length(traj) != length(cov)) {
        stop("If 'cov' is not 'NA' then its length (", length(cov), ") has to match that of 'traj' (", length(traj), ").")
      }
    }

    # only if the s-estimate is not 'NA' then try to estimate p-value for s-estimate
    if(!is.na(param["s"])) {

      # simulate AF trajectories with random drift only
      repl <- nrow(traj)
      simTraj <- wf.traj(p0=rep(p0, times=repl*N.pval), Ne=Ne, t=t, s=0, h=0.5, haploid=haploid)

      # add sampling noise according to 'cov'
      if(any(!is.na(cov))) {
        smpl <- sample.alleles(p=as.vector(t(simTraj)), size=as.vector(t(cov)), mode="coverage", ploidy=if(haploid) 1 else 2)
        simTraj <- matrix(smpl, ncol=ncol(simTraj), byrow=TRUE, dimnames=dimnames(simTraj))
      }

      # estimate s/h for each simulated locus
      sim.param <- foreach(l=1:N.pval, .combine=rbind) %do% {

        cT <- consensus.traj(traj=simTraj[seq((l-1)*repl+1, l*repl),], t=t, cov=cov, bias=cBias)

        if(method == "LLS" || (method == "automatic" && useLLS(ctraj=cT, h=h, haploid=haploid))) {
          # perofrm linear regression (approximation)
          lm.s(ctraj=cT, haploid=haploid)
        } else {
          # perform non-linear least-squares regression
          nls.sh(ctraj=cT, Ne=Ne, h=h, haploid=haploid, s.start=0.1, h.start=0.5, approximate=approximate)
        }
      }

      # compute p-value
      pval <- sum(abs(sim.param[,"s"]) >= abs(param["s"]), na.rm=TRUE) / sum(!is.na(sim.param[,"s"]))
      if(sum(is.na(sim.param[,"s"])) > N.pval*0.05)
        warning("P-value might be inaccurate under the given scenario.")
    }

    estsh.out$p.value <- pval
    estsh.out$N.pval <- N.pval
  }

  class(estsh.out) <- "estsh"
  return(estsh.out)
}

compareSH <- function(traj, t, Ne, haploid=FALSE, h=NA, N.ctraj=10000, N.pval=1000, cov=NA, method=c("LLS", "NLS", "automatic")) {
  if(is.vector(traj)) {
    traj <- matrix(traj, nrow=1)
  } else {
    traj <- as.matrix(traj)
  }
  if(ncol(traj) <= 1)
    stop("Allele frequencies at more than one time point have to be specified.")
  if(ncol(traj) != length(t <- as.numeric(t)))
    stop("Number of columns in 'traj' has to be equal to the length of 't'.")
  if(length(Ne <- as.numeric(Ne)) != 1)
    stop("Length of 'Ne' has to be equal to 1.")

  # common estimate of s/h for all replicates
  cparam <- estimateSH(traj=traj, t=t, Ne=Ne, haploid=haploid, h=h, N.ctraj=N.ctraj, simulate.p.value=FALSE, cov=cov, method=method)

  # separate estimate of s/h for each replicate
  sparam <- apply(traj, 1, function(trj) {
    coef(estimateSH(traj=trj, t=t, Ne=Ne, haploid=haploid, h=h, N.ctraj=N.ctraj, simulate.p.value=FALSE, cov=cov, method=method))
  })

  # if common estimate of s/h cannot be obtained than stop
  if(is.na(cparam$s)) {
    stop("Common estimate of s cannot be obtained.")
  }

  # generate result
  compsh.out <- list(traj=traj, t=t, Ne=Ne, haploid=haploid, N.ctraj=N.ctraj, N.pval=N.pval, s=cparam$s,
                     s.sep=if(is.vector(sparam)) sparam else sparam["s",])
  if(is.na(h)) {
    compsh.out$h <- cparam$h
    compsh.out$h.sep <- sparam["h",]
  } else {
    compsh.out$h.given <- h
  }

  # simulate AF trajectories assuming both s and h are identical for all replicates
  repl <- nrow(traj)
  p0 <- cparam$ctraj$af[1]
  simTraj <- wf.traj(p0=rep(p0, times=repl*N.pval), Ne=Ne, t=t, s=cparam$s, h=if(is.na(h)) cparam$h else h, haploid=haploid)

  # introduce sampling noise
  if(any(!is.na(cov))) {
    # make sure that 'cov' is set properly
    if(length(traj) != length(cov)) {
      stop("If 'cov' is not 'NA' then its length (", length(cov), ") has to match that of 'traj' (", length(traj), ").")
    }
    compsh.out$cov <- cov
    # add sampling noise according to 'cov'
    smpl <- sample.alleles(p=as.vector(t(simTraj)), size=as.vector(t(cov)), mode="coverage", ploidy=if(haploid) 1 else 2)
    simTraj <- matrix(smpl, ncol=ncol(simTraj), byrow=TRUE, dimnames=dimnames(simTraj))
  }

  # estimate s/h for each simulated locus separately
  sparam.sim <- apply(simTraj, 1, function(trj) {
    coef(estimateSH(traj=trj, t=t, Ne=Ne, haploid=haploid, h=h, N.ctraj=N.ctraj, simulate.p.value=FALSE, cov=cov))
  })

  # compute p-value for variance of s
  var.emp <- if(is.vector(sparam)) var(sparam, na.rm=TRUE) else var(sparam["s",], na.rm=TRUE)
  var.sim <- if(is.vector(sparam.sim)) colVars(matrix(sparam.sim, nrow=repl), na.rm=TRUE) else colVars(matrix(sparam.sim["s",], nrow=repl), na.rm=TRUE)
  compsh.out$s.p.value <- sum(var.sim >= var.emp, na.rm=TRUE) / sum(!is.na(var.sim))
  if(sum(is.na(compsh.out$s.p.value)) > N.pval*0.05)
    warning("P-value for 's' might be inaccurate under the given scenario.")

  # compute p-value for variance of h
  if(is.na(h)) {
    var.emp <- var(sparam["h",], na.rm=TRUE)
    var.sim <- colVars(matrix(sparam.sim["h",], nrow=repl), na.rm=TRUE)
    compsh.out$h.p.value <- sum(var.sim >= var.emp, na.rm=TRUE) / sum(!is.na(var.sim))
    if(sum(is.na(var.sim)) > N.pval*0.05)
      warning("P-value for 'h' might be inaccurate under the given scenario.")
  }

  class(compsh.out) <- "compsh"
  return(compsh.out)
}



