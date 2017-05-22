\name{nls.sh}
\alias{nls.sh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate s and h with NLS
}
\description{
\code{nls.sh} estimates s and h with non-linear least squares regression.
}
\usage{
nls.sh(ctraj, Ne, h, haploid, s.start = 0.1, h.start = 0.5, approximate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ctraj}{
a consensus trajectory as returned by \code{\link{consensus.traj}}.
}
  \item{Ne}{
numeric specifying the effective population size
}
  \item{h}{
numeric defining the dominance that should be assumed. If \code{h = NA} (default) dominance will be estimated alongside with the selection coefficient.
}
  \item{haploid}{
logical indicating wether a population of haploid (\code{haploid = TRUE}) or diploid (default) individuals should be simulated.
}
  \item{s.start}{
numeric starting value for s provided to \code{\link{nls}}.
}
  \item{h.start}{
numeric starting value for h provided to \code{\link{nls}}.
}
  \item{approximate}{
logical indicating if a continuous time approximation should be used (\code{approximate = TRUE}), or not (default).
}
}
\details{
Instead of calling this function directly, use \code{\link{estimateSH}} instead.
}
\value{
A named numeric vector with estimates of both s and h (if \code{h} equals \code{NA}) is returned.
}
\references{
Taus T., Futschik A. and Schlötterer C., manuscript in preparation.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
