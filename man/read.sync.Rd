\name{read.sync}
\alias{read.sync}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read Sync-File Format
}
\description{
Reads a synchronized file generated e.g. by the \href{http://sourceforge.net/projects/popoolation2/}{Popoolation2} software and extracts biallelic count data. Allele frequencies are polarized either by the minor or the rising allele.
}
\usage{
read.sync(file, gen, repl, rising = FALSE)
}
\arguments{
  \item{file}{
the name of the file, which the data is to be read from.
}
  \item{gen}{
numeric vector specifying the number of generations for each population, see 'Details'.
}
  \item{repl}{
numeric vector defining how populations should be grouped by replicate, see 'Details'.
}
  \item{rising}{
a logical value. If \code{FALSE} (default) allele counts are polarized for the minor allele. If \code{TRUE} allele counts are polarized for the rising allele.
}
}
\details{
MISSING
}
\value{
\code{read.sync} returns an object of class \code{sync}. Allelic information, allele frequencies and sequence coverages can be accessed with various functions including \code{\link{af.traj}}, \code{\link{af}} and \code{\link{coverage}}
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
\code{\link{af.traj}}, \code{\link{af}} and \code{\link{coverage}}
}
\examples{
#MISSING
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
