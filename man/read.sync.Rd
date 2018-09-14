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
read.sync(file, gen, repl, polarization = c("minor", "rising", "reference"), keepOnlyBiallelic = FALSE)
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
  \item{polarization}{
character determining, how allele frequencies should be polarized. If \code{polarization = "minor"} (default) allele counts are polarized for the minor allele. If \code{polarization = "rising"} or \code{polarization = "reference"} then allele counts are polarized for the rising or reference allele, respectively.
}
  \item{keepOnlyBiallelic}{
specifies if only polymorphic sites should be kept. Combining allele counts of all populations, monomorphic sites will always be removed. If \code{keepOnlyBiallelic = FALSE} (default), then also triallelic sites will be kept. Otherwise (\code{keepOnlyBiallelic = TRUE}), only biallelic sites will be kept.
}
}
\details{

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
\code{\link{alleles}}, \code{\link{polarization}}, \code{\link{af.traj}}, \code{\link{af}} and \code{\link{coverage}}
}
\examples{
#MISSING
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
