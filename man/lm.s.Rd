\name{lm.s}
\alias{lm.s}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate s with LLS
}
\description{
\code{lm.s} estimates s with linear least squares regression.
}
\usage{
lm.s(ctraj, haploid, maxiter = 10, tol = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ctraj}{
a consensus trajectory as returned by \code{\link{consensus.traj}}.
}
  \item{haploid}{
logical indicating wether a population of haploid (\code{haploid = TRUE}) or diploid (default) individuals should be simulated.
}
  \item{maxiter}{
integer defining the maximum number of iterations to correct for discrepancies between a selection model with discrete or overlapping generations.
}
  \item{tol}{
numeric specifying the tolerance for the iterative bias correction. If \code{tol = 0.01}, then the correction procedure will stop once the deviation between the discrete and continuouse model is less than 1\%.
}
}
\details{
Instead of calling this function directly, use \code{\link{estimateSH}} instead.
}
\value{
A named numeric vector with estimates of both s and p0 is returned.
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
\code{\link{estimateSH}} and \code{\link{nls.sh}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
