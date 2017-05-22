\name{splitLocusID}
\alias{splitLocusID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Split Locus ID
}
\description{
Splits the ID of a genomic locus into chromosome and position.
}
\usage{
splitLocusID(id, sep = ".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{
character vector providing locus IDs.
}
  \item{sep}{
a character string to separate chromosome name and position.
}
}
\details{

}
\value{
\code{splitLocusID} returns a \code{\link{data.table}} with the following columns
\item{chr}{chromosome name.}
\item{pos}{position.}
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
\code{\link{af}} and \code{\link{coverage}}.
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
