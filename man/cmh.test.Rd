\name{cmh.test}
\alias{cmh.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Vectorized Cochran-Mantel-Haenszel Test for Count Data
}
\description{
Performs a vectorized version of the Cochran-Mantel-Haenszel (CMH) test with the nullhypothesis that relative allele frequencies are not associated with time point.
}
\usage{
cmh.test(A0, a0, At, at, min.cov = 1, max.cov = 1, min.cnt = 1, log = FALSE, return.only.pval = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A0, a0, At, at}{
numeric matrix with biallelic count data at the first (\code{A0}, \code{a0}) and second (\code{At}, \code{at}) time point. Rows and columns correspond to replicates and individual loci, respectively.
}
  \item{min.cov}{
numeric specifying the minimal sequence coverage. If the sum of biallelic counts (\code{A}, \code{a}) are below \code{min.cov} in any of the replicates at any of the two time points (\code{0}, \code{t}), then \code{cmh.test} will return \code{NA} for the respective loci.
}
  \item{max.cov}{
numeric specifying the maximal sequence coverage. Values between \code{0} and \code{1} are interpreted as quantile thresholds, see 'Details'. If the sum of biallelic counts (\code{A}, \code{a}) exceeds the threshold in any of the replicates at any of the two time points (\code{0}, \code{t}), then \code{cmh.test} will return \code{NA} for the respective loci. By default (\code{max.cov = 1}) the CMH-test will be applied to all loci.
}
  \item{min.cnt}{
numeric indicating the minimal minor allele count. If the sum of minor allele counts over all replicates and time points is below \code{min.cnt}, then \code{cmh.test} will return \code{NA} for the respective loci.
}
  \item{log}{
  logical determining whether p-values should be returned directly (default), or after \code{-log10} transformation (if \code{log = TRUE}).
}
  \item{return.only.pval}{
  logical indicating whether to return only the p-values or p-values and test statistic in a named list.
  }
}
\details{
\code{max.cov} can be specified in two distinct ways. Values larger than 1 are interpreted as absolute sequence coverage thresholds, while values between 0 and 1 are interpreted as quantile thresholds. If for example \code{max.cov = 0.95} then \code{cmh.test} will return \code{NA} for all loci, where the sequence coverage exceeds the 95\% quantile in any of the replicates or time points. Please note that the quantile thresholds are computed separately for each replicate and time point.
}
\value{
\code{cmh.test} returns either a numeric vector of p-values (\code{return.only.pval = TRUE}),
or a named list of p-values and test statistic.
}
\references{
Agresti A.: Categorical data analysis (second edition). \emph{New York: Wiley} \bold{2002}
}
\author{
Thomas Taus
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{chi.sq.test}}
}
\examples{

# simulate data - 1 locus, 3 replicates
coverage <- 20
repl <- 3
p0 <- 0.1
pt <- 0.6
A0 <- rbinom(repl, size=coverage, prob=p0)
a0 <- coverage - A0
At <- rbinom(repl, size=coverage, prob=pt)
at <- coverage - At
# perform CMH-test
cmh.test(A0=A0, a0=a0, At=At, at=at)

# simulate data - 10 loci, 3 replicates
coverage <- 20
repl <- 3
loci <- 10
p0 <- 0.1
pt <- 0.6
A0 <- foreach(l=1:loci, .combine=cbind, .final=unname) \%do\% { rbinom(repl, size=coverage, prob=p0) }
a0 <- coverage - A0
At <- foreach(l=1:loci, .combine=cbind, .final=unname) \%do\% { rbinom(repl, size=coverage, prob=pt) }
at <- coverage - At
# show the structure of the simulated cound data
A0
a0
At
at
# perform CMH-test for each locus individually
cmh.test(A0=A0, a0=a0, At=At, at=at)

# get allele counts from empirical data
coverage0 <- t(coverage(dmelER, repl=1:3, gen=0))
A0 <- t(af(dmelER, repl=1:3, gen=0)) * coverage0
a0 <- coverage0 - A0
coverage59 <- t(coverage(dmelER, repl=1:3, gen=59))
A59 <- t(af(dmelER, repl=1:3, gen=59)) * coverage59
a59 <- coverage59 - A59
# perform CMH-test for all empirical loci
p.values <- cmh.test(A0=A0, a0=a0, At=A59, at=a59, min.cov=20, max.cov=0.99, min.cnt=1, log=TRUE)
# plot -log10 transformed p-values
plot(p.values, main="CMH Manhattan plot", xlab="SNP", ylab="-log10(p)", pch=".")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
