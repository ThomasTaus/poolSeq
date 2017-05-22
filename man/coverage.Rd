\name{coverage}
\alias{coverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sequence coverages
}
\description{
Extract sequence coverages of all genomic loci from a \code{sync} object.
}
\usage{
coverage(sync, repl, gen)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sync}{
a \code{sync} object.
}
  \item{repl}{
numeric vector indicating the replicate(s).
}
  \item{gen}{
numeric vector indicating the generation(s).
}
}
\details{
}
\value{
\code{coverage} returns a matrix of sequence coverages for the replicates \code{repl} and generations \code{gen} specified. Rows and columns correspond to genomic loci and replicate-generation pairs, respectively. Please note that the order of genomic loci in the sync-file is not preserved.
}
\references{
}
\author{
Thomas Taus
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{af}}, \code{\link{af.traj}} and \code{\link{splitLocusID}}.
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
