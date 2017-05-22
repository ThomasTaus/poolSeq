\name{consensus.traj}
\alias{consensus.traj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Consensus Allele Frequency Trajectory
}
\description{
Combine multiple allele frequency trajectories to get a consensus trajectory.
}
\usage{
consensus.traj(traj, t, cov, Ne, N.sim = 10000, bias = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traj}{
numeric matrix of allele frequency trajectories, where rows and columns correspond to individual replicates and time points, respectively.
}
  \item{t}{
numeric vector indicating the number of generations for each column in \code{traj}.
}
  \item{cov}{
numeric matrix of sequence coverages, where columns and rows correspond to measurement time points and replicates, respectively
}
  \item{Ne}{
numeric specifying the effective population size
}
  \item{N.sim}{
numeric determining the number of simulations that should be performed to correct for bias (by default \code{sim = 10000}).
}
  \item{bias}{
numeric specifying the bias correction, which is necessesary to account for conditioning on fixation. By default (\code{bias = NA}) \code{\link{fixationBias}} is called to compute the bias.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with class "ctraj" containing the following components:
\item{af}{consensus allele frequencies.}
\item{t}{time in generations.}
\item{cov}{combined sequence coverages (if available).}
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
p0 <- 0.01
s <- 0.1
h <- 0.5
Ne <- 1000
tp <- seq(0, 100, by=10)
alleleFreqs <- wf.traj(p0=rep(p0, times=100), Ne=Ne, t=tp, s=s, h=h)

# compute consensus trajectory
cTraj <- consensus.traj(alleleFreqs, t=tp, cov=NA, Ne=Ne)

# plot individual and consensus allele frequency trajectories, as well as the trajectory assuming an infinitely large population (no random drift)
plot(1, type="n", xlim=c(0, max(tp)), ylim=c(0, 1), main="Allele frequency trajectories with selection", xlab="Generation", ylab="Allele frequency")
for(r in 1:nrow(alleleFreqs)) {
  lines(tp, alleleFreqs[r,])
}
lines(tp, cTraj$af, col="red", lwd=2)
lines(tp, wf.traj(p0=p0, Ne=NA, t=tp, s=s, h=h), col="green", lwd=2)
legend("topleft", legend=c("replicate", "consensus", "deterministic"), col=c("black", "red", "green"), lwd=c(1, 2, 2), bty="n")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
