\name{compareSH}
\alias{compareSH}
\alias{print.compsh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare Selection Coefficients and Dominance between Replicates
}
\description{
\code{compareSH} tests the nullhypothesis that replicated allele frequency trajectories were generated with identical selection coefficients (s) and dominance (h).
}
\usage{
compareSH(traj, t, Ne, haploid = FALSE, h = NA, N.ctraj = 10000, N.pval = 1000, cov = NA, method = c("LLS", "NLS", "automatic"))

\method{print}{compsh}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object of class "compsh" as returned by \code{compareSH}.
}
  \item{traj}{
numeric matrix of allele frequency trajectories, where columns correspond to individual time points
}
  \item{t}{
numeric vector indicating the number of generations for each column in \code{traj}.
}
  \item{Ne}{
numeric specifying the effective population size
}
  \item{haploid}{
logical indicating wether a population of haploid (\code{haploid = TRUE}) or diploid (default) individuals should be simulated.
}
  \item{h}{
numeric defining the dominance that should be assumed. If \code{h = NA} (default) dominance will be estimated alongside with the selection coefficient.
}
  \item{N.ctraj}{
numeric determining the number of simulations that should be performed to correct for bias in the consensus trajectory (by default \code{sim = 10000}).
}
  \item{N.pval}{
integer; the number of simulations performed to estimate the p-value (by default \code{N.pval = 1000}).
}
  \item{cov}{
numeric matrix of sequence coverages, where columns and rows correspond to measurement time points and replicates, respectively
}
  \item{method}{
  character string defining the method that should be used to estimate s (and p0 or h).
}
  \item{\dots}{
currently not used.
}
}
\details{

}
\value{
A list with class "compsh" containing the following components:
\item{traj}{matrix with observed allele frequency trajectories.}
\item{t}{generations of each measurement time point.}
\item{Ne}{effective population size.}
\item{haploid}{logical for diploid or haploid individuals.}
\item{ctraj}{consensus trajectory used to estimate s and h.}
\item{N.ctraj}{number of simulations to correct 'ctraj' from bias.}
\item{s}{estimated common s.}
\item{s.sep}{individual estimates of s, one for each replicate. }
\item{h}{if \code{h == NA}, estimated common h.}
\item{h.sep}{individual estimates of h, one for each replicate.}
\item{h.given}{dominance defined in \code{h}.}
\item{cov}{matrix with observed sequence coverages.}
\item{s.p.value}{p-value for s.}
\item{h.p.value}{p-value for h.}
\item{N.pval}{number of simulations to estimate p-value(s).}
}

\author{
Thomas Taus
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{estimateSH}}
}
\examples{
# simulate allele frequency trajectories under constant selection coefficient (and dominance)
tp <- seq(0, 100, by=10)
Ne <- 1000
traj <- wf.traj(p0=rep(0.05, times=3), Ne=Ne, t=tp, s=0.1, h=0.5, haploid=FALSE)

# test the null hypothesis that selection coefficients are identical between replicates
compareSH(traj=traj, t=tp, Ne=Ne, haploid=FALSE, N.pval = 100, h=0.5) ##### takes about 15 sec #####
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
