\name{dmelER}
\alias{dmelER}
\docType{data}
\title{
Evolve and Resequencing Data
}
\description{
\code{dmelER} contains allele count data of an evolve and resequence (E&R) study conducted with \emph{Drosophila melanogaster}.
}
\usage{data("dmelER")}
\format{
  \code{dmelER} is a \code{sync} object.
}
\details{
Three replicated founder populations were set up based on 113 isofemale lines established from a natural Northern Portuguese population of \emph{D. melanogaster}. The three populations were propagated as independent replicates in a hot fluctuating laboratory environment, with 28°C and 18°C during light and dark periods, respectively. At multiple time points (F0, F15, F37 and F59) Pool-seq was applied to obtain genome-wide allele frequency estimates for each of the three replicates. Temporal allele frequencies are available for a subset of about 15,000 genomic loci.
}
\source{
Orozco-terWengel P., Kapun M., Nolte V., Kofler R., Flatt T., & Schlötterer C.: Adaptation of Drosophila to a novel laboratory environment reveals temporally heterogeneous trajectories of selected alleles, \emph{Molecular Ecology} \bold{2012}, 21, 4931–41.

Franssen S. U., Nolte V., Tobler R., & Schlötterer C.: Patterns of linkage disequilibrium and long range hitchhiking in evolving experimental Drosophila melanogaster populations, \emph{Molecular Biology and Evolution} \bold{2015}, 32, 495–509.
}
\seealso{
\code{\link{sync-class}}
}
\examples{
data(dmelER)
}
\keyword{datasets}
