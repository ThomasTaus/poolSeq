\name{wf.traj}
\alias{wf.traj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Random Genetic Drift
}
\description{
\code{wf.traj} simulates allele frequency trajectories based on the Wright-Fisher model.
}
\usage{
wf.traj(p0, Ne, t, s = 0, h = 0.5, haploid = FALSE, approximate = FALSE)
}
\arguments{
  \item{p0}{
numeric vector with starting allele frequencies.
}
  \item{Ne}{
integer specifying the desired effective population size. The unit is haploid (if \code{haploid = TRUE}) or diploid (if \code{haploid = FALSE}) individuals. If \code{is.na(Ne)} an infinitely large population size without random genetic drift will be simulated.
}
  \item{t}{
numeric vector defining for which generations allele frequencies should be returned.
}
  \item{s}{
  numeric specifying the selection coefficient used. By default (\code{s = 0}) no selection is applied.
}
  \item{h}{
  numeric indicating the dominance coefficient. This parameter will be ignored if \code{haploid = TRUE}.
}
  \item{haploid}{
  logical indicating wether a population of haploid (\code{haploid = TRUE}) or diploid (default) individuals should be simulated.
}
  \item{approximate}{
  logical. If \code{approximate = true} a continuous time approximation is used if applicable. Alternatively, forward in time simulations are performed (default).
}
}
\value{
\code{wf.traj} returns either a vector (if \code{length(p0) == 1} and \code{length(s) == 1} and \code{length(h) == 1}) or a matrix that contains allele frequency trajectories. Columns correspond to the time points defined in \code{t} and rows to individual replicates. If only one replicate should be simulated a numeric vector is returned instead of a matrix.

Allele frequencies in the output are guaranteed to be ordered increasing by the number of generations.
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
# simulate allele frequency trajectories of individual loci using the Wright-Fisher model (diploid individuals)
Ne <- 200
gen <- seq(0, 100, by=10)
alleleFreqs <- wf.traj(p0=rep(0.5, times=500), Ne=Ne, t=gen)
# look at a subset of the generated data
head(alleleFreqs)

# plot allele frequency trajectories
plot(1, type="n", xlim=c(0, max(gen)), ylim=c(0, 1), main="Neutral Genetic Drift", xlab="Generation", ylab="Allele frequency (\%)")
for(r in 1:nrow(alleleFreqs)) {
  lines(gen, alleleFreqs[r,])
}

# simulate allele frequency trajectories including selection
alleleFreqs <- wf.traj(p0=rep(0.05, times=500), Ne=Ne, t=gen, s=0.1, h=0.5)

# plot results
plot(1, type="n", xlim=c(0, max(gen)), ylim=c(0, 1), main="Positive Selection", xlab="Generation", ylab="Allele frequency (\%)")
for(r in 1:nrow(alleleFreqs)) {
  lines(gen, alleleFreqs[r,])
}

# add the trajectory under selection for a population of infinite size (no random genetic drift)
lines(gen, wf.traj(p0=0.05, Ne=NA, t=gen, s=0.1, h=0.5), col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
