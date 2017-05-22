\name{alleles}
\alias{alleles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Allelic information
}
\description{
Extract allelic information, such as reference, minor, major and rising allele, of all genomic loci from a \code{sync} object.
}
\usage{
alleles(sync)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sync}{
a \code{sync} object.
}
}
\details{
}
\value{
\code{alleles} returns a \code{\link{data.table}} containing information for each genomic locus in \code{sync}. Please note that the order of genomic loci in the sync-file is not preserved. The returned \code{data.table} provides the following columns:
\item{chr}{chromosome name.}
\item{pos}{sequence position.}
\item{ref}{reference allele.}
\item{major}{major allele.}
\item{minor}{minor allele.}
\item{rising}{rising allele.}
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
\code{\link{af}} and \code{\link{af.traj}}.
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
