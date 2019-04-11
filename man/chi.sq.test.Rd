\name{chi.sq.test}
\alias{chi.sq.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Vectorized Pearson's Chi-squared Test for Count Data
}
\description{
Performs vectorized tests of independence and goodness-of-fit tests for allele count data.
}
\usage{
chi.sq.test(A0, a0, At = NULL, at = NULL, p0 = 0.5, min.cov = 1, max.cov = 1, min.cnt = 1, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A0, a0, At, at}{
numeric vector with biallelic count data at the first (\code{A0}, \code{a0}) and second (\code{At}, \code{at}) time point. The i-th element of each vector is treated as count data of an individual locus.
}
  \item{p0}{
a numeric vector with expected relative frequencies of \code{A0}. Ignored if both \code{At} and \code{at} are given, see 'Details'.
}
 \item{min.cov}{
numeric specifying the minimal sequence coverage. If the sum of biallelic counts (\code{A}, \code{a}) are below \code{min.cov} in any of the replicates at any of the two time points (\code{0}, \code{t}), then \code{chi.sq.test} will return \code{NA} for the respective loci.
}
  \item{max.cov}{
numeric specifying the maximal sequence coverage. Values between \code{0} and \code{1} are interpreted as quantile thresholds, see 'Details'. If the sum of biallelic counts (\code{A}, \code{a}) exceeds the threshold in any of the replicates at any of the two time points (\code{0}, \code{t}), then \code{chi.sq.test} will return \code{NA} for the respective loci. By default (\code{max.cov = 1}) the chi-squared test will be applied to all loci.
}
  \item{min.cnt}{
numeric indicating the minimal minor allele count. If the sum of minor allele counts over all replicates and time points is below \code{min.cnt}, then \code{chi.sq.test} will return \code{NA} for the respective loci.
}
  \item{log}{
  logical determining whether p-values should be returned directly (default), or after \code{-log10} transformation (if \code{log = TRUE}).
}
  \item{return.only.pval}{
  logical indicating whether to return only the p-values or p-values and test statistic in a named list.
  }
}
\details{
If either \code{At} or \code{at} is missing, then a \emph{goodness-of-fit test} is performed with the nullhypothesis that relative frequencies of \code{A}-alleles are equal to those specified in \code{p0}. If both \code{At} and \code{at} are given, then \emph{Pearson's chi-squared test of independence} is performed with the nullhypothesis that relative allele frequencies are not associated with time point.

\code{max.cov} can be specified in two distinct ways. Values larger than 1 are interpreted as absolute sequence coverage thresholds, while values between 0 and 1 are interpreted as quantile thresholds. If for example \code{max.cov = 0.95} then \code{chi.sq.test} will return \code{NA} for all loci, where the sequence coverage exceeds the 95\% quantile in any of the replicates or time points. Please note that the quantile thresholds are computed separately for each replicate and time point.
}
\value{
\code{chi.sq.test} returns either a numeric vector of p-values (\code{return.only.pval = TRUE}),
or a named list of p-values and test statistic.
}
\references{
Sokal R. S. and Rohlf F. J.: Biometry (third edition). \emph{New York: W. H. Freeman and Company} \bold{1995}
}
\author{
Thomas Taus
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cmh.test}}
}
\examples{
# simulate allele counts - 2 time points, 10 individual loci
coverage <- 20
loci <- 10
p0 <- 0.1
pt <- 0.6
A0 <- rbinom(loci, size=coverage, prob=p0)
a0 <- coverage - A0
At <- rbinom(loci, size=coverage, prob=pt)
at <- coverage - At
# perform Pearson's chi-squared test of independence
chi.sq.test(A0=A0, a0=a0, At=At, at=at)

# simulate allele counts - 1 time point, 10 individual loci
coverage <- 20
loci <- 10
p0 <- 0.2
A0 <- rbinom(loci, size=coverage, prob=p0)
a0 <- coverage - A0
# perform goodness-of-fit test
chi.sq.test(A0=A0, a0=a0, p0=0.5)

parOld <- par(mfrow=c(3,1))
for(r in 1:3) {
  # get allele counts from empirical data
  coverage0 <- t(coverage(dmelER, repl=r, gen=0))
  A0 <- t(af(dmelER, repl=r, gen=0)) * coverage0
  a0 <- coverage0 - A0
  coverage59 <- t(coverage(dmelER, repl=r, gen=59))
  A59 <- t(af(dmelER, repl=r, gen=59)) * coverage59
  a59 <- coverage59 - A59
  # apply Pearson's chi-squared test of independence to each replicate separately
  p.values <- chi.sq.test(A0=A0, a0=a0, At=A59, at=a59, min.cov=20, max.cov=0.99, min.cnt=1, log=TRUE)
  # plot -log10 transformed p-values
  plot(p.values, main=paste0("Chi-squared test - Replicate ", r), ylim=c(0, max(p.values, na.rm=TRUE)), xlab="SNP", ylab="-log10(p)", pch=".")
}
par(parOld)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
