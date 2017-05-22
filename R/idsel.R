# -----------------------------------------------------
# Statistical tests to contrast allele frequency data -
# -----------------------------------------------------

cmh.test <- function(A0, a0, At, at, min.cov = 1, max.cov = 1, min.cnt = 1, log = FALSE) {
  # make sure that dimensions of A0, a0, At and at are identical
  if(!identical(dim(A0 <- as.matrix(A0)), dim(a0 <- as.matrix(a0))) ||
     !identical(dim(A0), dim(At <- as.matrix(At))) ||
     !identical(dim(A0), dim(at <- as.matrix(at)))) {
    stop("Dimension of 'A0' (", dim(A0), "), 'a0' (", dim(a0), "), 'At' (", dim(At), ") and 'at' (", dim(at), ") are not identical.")
  }
  # make sure that there are at least 2 replicates
  if(nrow(A0) < 2) {
    stop("You have to provide count data for at least 2 replicates.")
  }
  # make sure that A0, a0, At and at are numeric
  if(!is.numeric(A0) | !is.numeric(a0) | !is.numeric(At) | !is.numeric(at)) {
    stop("Either 'A0', 'a0', 'At' or 'at' is not numeric")
  }
  # make sure there are no negative values
  if(any(A0 < 0 | a0 < 0 | At < 0 | at < 0, na.rm=TRUE)) {
    stop("Negative counts are not allowed in neither 'A0', 'a0', 'At' nor 'at'.")
  }
  # make sure that min.cov, max.cov, min.cnt and log are set properly
  if(length(min.cov <- as.numeric(min.cov)) != 1 | length(max.cov <- as.numeric(max.cov)) != 1 | length(min.cnt <- as.numeric(min.cnt)) != 1 | length(log <- as.logical(log)) != 1) {
    stop("Length of 'min.cov' (", length(min.cov), "), 'max.cov' (", length(max.cov), "), 'min.cnt' (", length(min.cnt), ") and 'log' (", length(log), ") has to be equal to '1'.")
  }
  if(is.na(min.cov) | min.cov < 1 | is.na(min.cnt) | min.cnt < 1) {
    stop("Both 'min.cov' (", min.cov, ") and 'min.cnt' (", min.cnt, ") have to be >= 1.")
  }
  if(is.na(max.cov) | max.cov < 0) {
    stop("'max.cov' (", max.cov, ") cannot be 'NA' or negative.")
  }
  # mask loci that do not meet either 'min.cov', 'max.cov' or 'min.cnt' requirements
  cov0 <- A0+a0
  covt <- At+at
  max.cov0 <- if(max.cov <= 1) apply(cov0, 1, quantile, probs=max.cov, na.rm=TRUE) else max.cov
  max.covt <- if(max.cov <= 1) apply(covt, 1, quantile, probs=max.cov, na.rm=TRUE) else max.cov
  mask <- colAnys(cov0 < min.cov) | colAnys(covt < min.cov) | colAnys(cov0 > max.cov0) | colAnys(covt > max.covt) | colSums(A0+At) < min.cnt | colSums(a0+at) < min.cnt
  A0[,mask] <- NA

  # compute CMH statistic
  n <- A0+a0+At+at
  CMH.chi <- (abs(colSums(A0-(A0+a0)*(A0+At)/n))-0.5)^2 / (colSums((A0+a0)*(A0+At)*(a0+at)*(At+at)/(n^3-n^2)))

  # return p-values according to chi-squared distribution
  p.value <- pchisq(unname(CMH.chi), df=1, lower.tail=FALSE)
  return(if(log) -log10(p.value) else p.value)
}

chi.sq.test <- function(A0, a0, At = NULL, at = NULL, p0 = 0.5, min.cov = 1, max.cov = 1, min.cnt = 1, log = FALSE) {
  # make sure that lengths of A0 and a0 are identical
  if(!identical(length(A0 <- as.numeric(A0)), length(a0 <- as.numeric(a0)))) {
    stop("Lengths of 'A0' (", length(A0), ") and 'a0' (", length(a0), ") are not identical.")
  }
  # make sure there are no negative values
  if(any(A0 < 0 | a0 < 0, na.rm=TRUE)) {
    stop("Negative counts are not allowed in neither 'A0' nor 'a0'.")
  }
  # make sure that min.cov, max.cov, min.cnt and log are set properly
  if(length(min.cov <- as.numeric(min.cov)) != 1 | length(max.cov <- as.numeric(max.cov)) != 1 | length(min.cnt <- as.numeric(min.cnt)) != 1 | length(log <- as.logical(log)) != 1) {
    stop("Length of 'min.cov' (", length(min.cov), "), 'max.cov' (", length(max.cov), "), 'min.cnt' (", length(min.cnt), ") and 'log' (", length(log), ") has to be equal to '1'.")
  }
  if(is.na(min.cov) | min.cov < 1 | is.na(min.cnt) | min.cnt < 1) {
    stop("Both 'min.cov' (", min.cov, ") and 'min.cnt' (", min.cnt, ") have to be >= 1.")
  }
  if(is.na(max.cov) | max.cov < 0) {
    stop("'max.cov' (", max.cov, ") cannot be 'NA' or negative.")
  }
  # if only one of either 'At' or 'at' is provided then throw a warning
  if((!is.null(At) && is.null(at)) || (is.null(At) && !is.null(at))) {
    warning("Only one of either 'At' or 'at' is provided. A goodness-of-fit test will be performed.")
  }

  # if both 'At' and 'at' are available then perform a TEST OF INDEPENDENCE
  if (!is.null(At) && !is.null(at)) {
    # make sure that lengths of A0, a0 ,At and at are identical
    if(!identical(length(A0), length(At <- as.numeric(At))) || !identical(length(A0), length(at <- as.numeric(at)))) {
      stop("Lengths of 'A0' (", length(A0), "), 'a0' (", length(a0), "), 'At' (", length(At), ") and 'at' (", length(at), ") are not identical.")
    }
    # make sure there are no negative values
    if(any(At < 0 | at < 0, na.rm=TRUE)) {
      stop("Negative counts are not allowed in neither 'A0' nor 'a0'.")
    }
    # mask loci that do not meet either 'min.cov' or 'min.cnt' requirements
    cov0 <- A0+a0
    covt <- At+at
    max.cov0 <- if(max.cov <= 1) quantile(cov0, probs=max.cov, na.rm=TRUE) else max.cov
    max.covt <- if(max.cov <= 1) quantile(covt, probs=max.cov, na.rm=TRUE) else max.cov
    mask <- cov0 < min.cov | covt < min.cov | cov0 > max.cov0 | covt > max.covt | A0+At < min.cnt | a0+at < min.cnt
    A0[mask] <- NA
    # compute test statistic
    n <- A0+a0+At+at
    x <- (abs(A0*at-a0*At)-n/2)^2*n / ((A0+a0)*(At+at)*(A0+At)*(a0+at))
    # return p-values according to chi-squared distribution
    p.value <- pchisq(x, df=1, lower.tail=FALSE)
    return(if(log) -log10(p.value) else p.value)

    # otherwise perform a GOODNESS-OF-FIT TEST (assuming 'A0' to occur at a relative frequency of 'p0')
  } else {
    # make sure that 'p0' is between 0 and 1
    if(any(is.na(p0) | p0 < 0 | p0 > 1)) {
      stop("'p0' has to be within the range [0, 1].")
    }
    # recycle p0
    if(length(A0) %% length(p0) != 0) {
      stop("Length of 'p0' (", length(p0), ") is not a multiple of the length of 'A0' (", length(A0), ").")
    }
    p0 <- rep_len(p0, length(A0))
    # mask loci that do not meet either 'min.cov' or 'min.cnt' requirements
    cov0 <- A0+a0
    max.cov0 <- if(max.cov <= 1) quantile(cov0, probs=max.cov, na.rm=TRUE) else max.cov
    mask <- cov0 < min.cov | cov0 > max.cov0 | A0 < min.cnt | a0 < min.cnt
    A0[mask] <- NA
    # compute test statistic
    n <- A0+a0
    x <- (A0-p0*n)^2 / (p0*(1-p0)*n)
    # return p-values according to chi-squared distribution
    p.value <- pchisq(x, df=1, lower.tail=FALSE)
    return(if(log) -log10(p.value) else p.value)
  }
}
