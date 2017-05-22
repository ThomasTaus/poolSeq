\name{af}
\alias{af}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Allele frequencies
}
\description{
Extract allele frequencies of all genomic loci from a \code{sync} object.
}
\usage{
af(sync, repl, gen)
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
\code{af} returns a matrix of allele frequencies for the replicates \code{repl} and generations \code{gen} specified. Rows and columns correspond to genomic loci and replicate-generation pairs, respectively. Please note that the order of genomic loci in the sync-file is not preserved.
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
\code{\link{af.traj}}, \code{\link{coverage}} and \code{\link{splitLocusID}}.
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
