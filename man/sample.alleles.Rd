\name{sample.alleles}
\alias{sample.alleles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample alleles
}
\description{
Introduce sampling variance in allele frequency data to mimic the Pool-seq approach. On one hand, subjecting only a subset of individuals in a population to Pool-seq is modeled with hypergeometric sampling (\code{mode = "individuals"}). On the other hand, sampling variance introduced by sequencing only a fraction of all DNA fragments is modeled with binomial sampling (\code{mode = "coverage"}).
}
\usage{
sample.alleles(p, size, mode = c("coverage", "individuals"), Ncensus = NA, ploidy = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
numeric vector defining relative allele frequencies, which are used as success probabilities in the sampling process.
  }
  \item{size}{
numeric indicating the sample size to be used for binomial (\code{mode = "coverage"}) or Poisson sampling (\code{mode = "individuals"}), see 'Details'.
  }
  \item{mode}{
character string specifying the sampling mode. Possible values are \code{"coverage"} and \code{"individuals"}.
  }
  \item{Ncensus}{
numeric specifying the census size of the entire population (before sampling).
  }
  \item{ploidy}{
numeric, the ploidy of the individuals.
  }
}
\details{
If \code{mode = "coverage"} and \code{length(size) == 1} then for each allele frequency an individual sequence coverage value will be drawn from a Poisson distribution with \code{lambda = size}. Otherwise (\code{length(size) > 1}) the values in \code{size} will be used directly and recycled if necessary. The \code{"coverage"} sampling mode applies \code{rbinom} with \code{size} equal to the sequence coverage and \code{prob} equal to the allele frequency (\code{p}).

If \code{mode = "individuals"} then \code{size} has to be an integer specifying the number of individuals with a certain \code{ploidy} that are sampled from the population. Here \code{rhyper} is applied.
}
\value{
A numeric vector of allele frequencies after introducing sampling variance or (if \code{mode = "coverage"} and \code{length(size) == 1}) a \code{data.table} containing the following columns: 
  \item{p.smpld}{allele frequencies after sampling}
  \item{size}{sequence coverage for each position, drawn from a Poisson distribution with \code{lambda = size}}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Thomas Taus
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# generate random allele frequencies
af <- runif(10000, min=0, max=1)

# introduce sampling variance to mimic Pool-seq of the entire population at 100X coverage
afSeq <- sample.alleles(af, size=100, mode="coverage")

# plot distribution of differences in allele frequency before and after sampling
hist(af-afSeq$p.smpld, main="Sequencing at 100X", xlab="Error in allele frequency (\%)", ylab="Occurrences")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
