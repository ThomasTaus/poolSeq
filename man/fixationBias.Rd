\name{fixationBias}
\alias{fixationBias}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bias Correction
}
\description{
\code{fixationBias} computes the bias that occurs, when conditioning on fixation of an allele.
}
\usage{
fixationBias(traj, t, Ne, N.sim = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traj}{
numeric matrix of allele frequency trajectories, where rows and columns correspond to individual replicates and time points, respectively.
}
  \item{t}{
numeric vector indicating the number of generations for each column in \code{traj}.
}
  \item{Ne}{
numeric specifying the effective population size
}
  \item{N.sim}{
numeric determining the number of simulations that should be performed to correct for bias (by default \code{sim = 10000}).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Numeric vector with expected bias at each time point.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
