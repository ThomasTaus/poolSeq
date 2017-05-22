\name{estimateSH}
\alias{estimateSH}
\alias{print.estsh}
\alias{coef.estsh}
\alias{confint.estsh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Selection Coefficient and Dominance
}
\description{
Based on a set of allele frequency trajectories and the effective population size both the selection coefficient (s) and dominance (h) can be estimated.
}
\usage{
estimateSH(traj, t, Ne, haploid = FALSE, h = NA, N.ctraj = 0, simulate.p.value = FALSE, N.pval = 1000, cov = NA, approximate = TRUE, method=c("LLS", "NLS", "automatic"))

\method{print}{estsh}(x, \dots)

\method{coef}{estsh}(object, \dots)

\method{confint}{estsh}(object, parm, level = 0.95, N.ci = 1000, warn = FALSE, \dots)
}
\arguments{
  \item{x, object}{
object of class "estsh" as returned by \code{estimateSH}.
}
  \item{traj}{
numeric matrix of allele frequency trajectories, where columns and rows correspond to measurement time points and replicates, respectively.
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
numeric determining the number of simulations that should be performed to correct for bias in the consensus trajectory if both p0 and Ne are small. By default, no bias correction is performed (\code{sim = 0}). If a bias correction should be performed, a value of \code{sim = 1000} is recommended.
}
  \item{simulate.p.value}{
locical indicating whether the p-value should be estimated by simulations.
}
  \item{N.pval}{
integer, the number of simulations performed to estimate the p-value (by default \code{N.pval = 1000}).
}
  \item{cov}{
numeric matrix of sequence coverages, where columns and rows correspond to measurement time points and replicates, respectively
}
  \item{approximate}{
  logical. If \code{approximate = true} a continuous time approximation is used to estimate s (and h). Alternatively, forward in time simulations are performed. The differences between continuous and discrete models increase with s.
}
  \item{method}{
  character string defining the method that should be used to estimate s (and p0 or h).
}
  \item{parm}{
a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
}
  \item{level}{
the confidence level required
}
  \item{N.ci}{
integer specifying the number of simulations that should be performed to get the bootstrap confidence interval.
}
  \item{warn}{
logical specifying if a warning should be sent in case the confidence interval may be inaccurate. Otherwise (default) \code{NA} will be returned in such cases.
}
  \item{\dots}{
currently not used.
}
}
\details{

}
\value{
A list with class "estsh" containing the following components:
\item{traj}{matrix with observed allele frequency trajectories.}
\item{t}{generations of each measurement time point.}
\item{Ne}{effective population size.}
\item{haploid}{logical for diploid or haploid individuals.}
\item{ctraj}{consensus trajectory used to estimate s and h.}
\item{N.ctraj}{number of simulations to correct 'ctraj' from bias.}
\item{s}{estimated selection coefficient.}
\item{h}{if \code{h == NA}, estimated dominance.}
\item{h.given}{dominance defined in \code{h}.}
\item{cov}{matrix with observed sequence coverages.}
\item{p.value}{p-value estimated by simulations under neutrality.}
\item{N.pval}{number of simulations to estimate the p-value.}
}
\author{
Thomas Taus
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{consensus.traj}} and \code{\link{compareSH}}
}
\examples{
# simulate allele frequency trajectories for a given selection coefficient and dominance
tp <- seq(0, 100, by=10)
Ne <- 1000
traj <- wf.traj(p0=rep(0.05, times=10), Ne=Ne, t=tp, s=0.1, h=0.5, haploid=FALSE)

# estimate only selection coefficient, assuming that dominance is 0.5
estimateSH(traj, t=tp, Ne=Ne, haploid=FALSE, h=0.5)

# estimate both selection coefficient and dominance for simulated data
estimateSH(traj, t=tp, Ne=Ne, haploid=FALSE, method="NLS")

# compute p-value for estimated selection coefficient
estimateSH(traj, t=tp, Ne=Ne, haploid=FALSE, h=0.5, simulate.p.value=TRUE) ##### takes about 15 sec #####

# estimate confidence interval for both s and h
res <- estimateSH(traj, t=tp, Ne=Ne, haploid=FALSE, h=0.5)
confint(res)  ##### takes about 15 sec #####
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
