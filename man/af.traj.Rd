\name{af.traj}
\alias{af.traj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Allele frequency trajectories
}
\description{
Extract allele frequency trajectories of specific genomic loci from a \code{sync} object.
}
\usage{
af.traj(sync, chr, pos, repl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sync}{
a \code{sync} object.
}
  \item{chr}{
character vector of chromosome names.
}
  \item{pos}{
numeric vector of sequence positions.
}
  \item{repl}{
numeric vector indicating the replicate(s).
}
}
\details{
}
\value{
\code{af.traj} returns either a matrix or a list of matrices containing allele frequencies for the loci (\code{chr} and \code{pos}) and replicates \code{repl} specified. Rows and columns correspond to genomic loci in a certain replicate and time points, respectively. 
}
\references{
}
\author{
Thomas Taus
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{af}} and \code{\link{coverage}}.
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
