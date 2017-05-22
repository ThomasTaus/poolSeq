\name{checkSNP}
\alias{checkSNP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtaine Appropriate SNPs for Ne Estimation
}
\description{
\code{checkSNP} helps to filter for those SNPs that are best suited for estimating Ne.
}
\usage{
checkSNP(p0, pt, cov0, covt, truncAF = NA)
}
\arguments{
  \item{p0}{
numeric vector of relative starting allele frequencies (\code{0 <= p0 <= 1}).
}
  \item{pt}{
numeric vector of relative allele frequencies after \code{t} generations (\code{0 <= pt <= 1}).
}
  \item{cov0}{
numeric vector of sequencing coverage at the first time point.
}
  \item{covt}{
numeric vector of sequencing coverage at the second time point.
}
  \item{truncAF}{
numeric indicating whether SNPs with extreme allele frequencies should be excluded for Ne estimation. If \code{truncAF = NA} (default) even extreme allele frequencies will be considered. If \code{truncAF = F} then only SNPs with an allele frequency \code{>= F} and \code{<= 1-F} will be considered for Ne estimation.
}
}
\details{
Alleles that start from a frequency of \code{0} or \code{1} are not informative for Ne estimation and \code{checkSNP} will return \code{FALSE} for those ones. Similarly positions with a sequencing coverage of 0 will be masked, and \code{checkSNP} will return \code{FALSE}.

Extreme allele frequencies, which are close to \code{0} or \code{1}, increase the variance of Ne estimates and it might therefore be advantageous not to consider them for Ne estimation. The parameter \code{truncAF} can be used for this purpose.
}
\value{
A logical vector indicating which SNPs are suited for Ne estimation.
}
\author{
Thomas Taus
}
\seealso{
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
