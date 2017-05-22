#-----------------------------------------
# Simulate allele frequency trajectories -
#-----------------------------------------

wf.traj <- function(p0, Ne, t, s=0, h=0.5, haploid=FALSE, approximate=FALSE) {

  max.len <- max(length(p0), length(s), length(h))
  traj <- NULL

  # if simulations should be performed for infinite size populations and continuous-time approximations are available
  if(is.na(Ne) && approximate && (haploid || all(h == 0.5))) {

    t.all <- rep(t, each=max.len)
    traj <- matrix(1/(1+(1-p0)/p0*exp(-s*t.all/if(haploid) 1 else 2)), ncol=length(t), dimnames=list(c(), paste0("F", t)))

    # perform forward in time simulations
  } else {

    # initialize trajectory matrix
    traj <- matrix(NA, ncol=length(t), nrow=max.len, dimnames=list(c(), paste0("F", t)))
    if(0 %in% t)
      traj[,"F0"] <- p0

    if(!haploid)
      Ne <- 2*Ne

    g <- 1
    p <- p0
    q <- 1-p0
    wAA <- 1+s
    wAa <- 1+h*s
    waa <- 1

    # simulate allele frequencies across time
    while(g <= max(t)) {
      # compute mean fitness
      w <- if(haploid) p*wAA + q*waa else p^2*wAA + 2*p*q*wAa + q^2*waa

      # apply selection and random drift
      p <- if(haploid) p*wAA/w else (wAA*p^2 + wAa*p*q)/w
      if(!is.na(Ne))
        p <- rbinom(length(p), Ne, p) / Ne
      q <- 1-p

      # if necessary then save current allele frequency to results
      if(g %in% t)
        traj[,paste0("F", g)] <- p

      g <- g+1
    }
  }

  if(nrow(traj) == 1)
    return(as.vector(traj))

  return(traj)
}

sample.alleles <- function(p, size, mode=c("coverage", "individuals"), Ncensus=NA, ploidy=2) {
  # determin number of return values
  maxlen <- max(length(p), length(size))
  # check length of a, n and size parameter
  if(maxlen %% length(p) != 0 || maxlen %% length(size) != 0)
    warning("Parameters differ in length and are not multiple of one another: length(p)=", length(p), ", length(size)=", length(size))

  # check mode parameter
  mode <- match.arg(mode)

  # sample allele counts according to the specified method
  return(switch(mode,
                coverage={
                  # if length of 'size' equals '1' then generate target coverage values using the Poisson distribution, otherwise use values of 'size' directly
                  cov <- if(length(size) == 1) rpois(n=maxlen, lambda=size) else size
                  # sample allele frequencies from Binomial distribution based on 'cov' and 'p'
                  p.smpld <- rbinom(n=maxlen, size=cov, prob=p) / cov
                  # return results, including coverage values if they were drawn from a Poisson distribution
                  if(length(size) == 1) data.table(p.smpld=p.smpld, size=cov) else p.smpld
                },
                individuals={
                  # if length of 'size' is larger than 1, then send warning message
                  if(length(size) > 1)
                    warning("Only the first element in 'size' will be used, because sampling mode is 'individuals'")

                  # sample random allele frequencies from Hypergeometric distribution
                  rhyper(nn=maxlen, m=p*Ncensus*ploidy, n=(1-p)*Ncensus*ploidy, k=size[1]*ploidy)/(size[1]*ploidy)
                }))
}

