\name{estimateWndNe}
\alias{estimateWndNe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Effective Population Size for Genome-Wide Data
}
\description{
\code{estimateNe} splits the genome-wide data into non-overlapping windows of size \code{wndSize} and estimates for each window the effective population size (Ne) based on temporal allele frequency changes.
}
\usage{
estimateWndNe(chr, pos, wndSize, p0, pt, cov0, covt, t, unit = c("bp", "SNP"), ploidy = 2, truncAF = NA, method = "P.planI", Ncensus = NA, poolSize = rep(Ncensus, times = 2))
}
\arguments{
  \item{chr}{
vector of chromosome identifiers.
}
  \item{pos}{
numeric vector of sequence positions.
}
  \item{wndSize}{
integer indicating the window size to be used for Ne estimation. If \code{NA}, Ne will be estimated for entire chromosomes. The interpretation of this parameter dependes on how \code{unit} is specified.
}
  \item{p0}{
numeric vector of relative starting allele frequencies (\code{0 <= p0 <= 1}).
}
  \item{pt}{
numeric vector of relative allele frequencies after \code{t} generations (\code{0 <= pt <= 1}).
}
  \item{cov0}{
numeric vector of sequencing coverage at the first time point.
}
  \item{covt}{
numeric vector of sequencing coverage at the second time point.
}
  \item{t}{
integer specifying the number of generations between the two allele frequency measurements.
}
  \item{unit}{
character string indicating the unit of \code{wndSize}. Possible values are \code{"bp"} (number of base pairs) and \code{"SNP"} (number of single-nucleotide polymorphisms).
}
  \item{ploidy}{
integer specifying the ploidy of the organisms in the population. Diploids (\code{ploidy = 2}) are the default.
}
  \item{truncAF}{
numeric indicating whether SNPs with extreme allele frequencies should be excluded from Ne estimation, see \code{\link{checkSNP}}.
}
  \item{method}{
character vector defining the method(s) that should be used to estimate Ne. Possible values include \code{"P.planI"}, \code{"P.planII"}, \code{"JR.planI"}, \code{"JR.planII"}, \code{"W.planI"} and \code{"W.planII"}, see 'Details'.
}
  \item{Ncensus}{
integer specifying the census population size, which is required only if \code{method} contains planI estimators.
}
 \item{poolSize}{
numeric vector of length 2 indicating the number of individuals pooled at the first and second time point. If the entire population has been subjected to Pool-seq (default) then \code{poolSize = rep(Ncensus, times=2)}, see 'Details'.
}
}
\details{
Different moment-based methods to estimate short-term Ne based on allele frequency changes between temporal samples are available. All estimators, whose names are passed on to \code{method} will be applied.
Plan I methods assume that individuals are sampled after reproduction or returned to the population after determining their genotypes, which results in a correlation of allele frequencies between the two sampling time points. An estimate of the census population size \code{Ncensus} is required.

In contrast, under plan II scenarios individuals are sampled before reproduction and their genotypes do not contribute to the next generation. Alleles frequencies at the two sampling time points will be uncorrelated and can be treated as two independent binomial samples.

Based on these two sampling schemes for individuals, estimators suggested by Waples 1989 (\code{"W.planI"} and \code{"W.planII"}), Jorde & Ryman 2007 (\code{"JR.planI"} and \code{"JR.planII"}) and Jónás et al. (\code{"P.planI"} and \code{"P.planII"}) are implemented. The latter are specifically design for application to Pool-seq data and, unlike the others, account for the two stage sampling process associated with such data.

The two sampling steps invovled in Pool-seq are: (i) individuals are drawn from the population without replacement (\code{poolSize}); (ii) reads are sampled out of the DNA pool (\code{cov0} and \code{covt}). Please note that all two-step sampling estimators require \code{poolSize} to be specified and unequal to \code{NA}.

In addition to the estimators mentioned above, Andreas Futschik suggested alternative estimators for one-step sampling (\code{"P.alt.1step.planII"}) and two-step sampling (\code{"P.alt.2step.planI"} and \code{"P.alt.2step.planII"}) that can be applied to Pool-seq data.
}
\value{
\code{estimateWndNe} returns an object of type \code{data.table} that contains Ne estimates for all windows of each chromosome using the method(s) specified in \code{method}. The resulting \code{data.table} contains the following columns:
  \item{Nx}{Ne estimate for the current window using method x, see 'Details'.}
  \item{chr}{chromosome name.}
  \item{wndStart}{}
  \item{wndStop}{start/stop position of the current window.}
  \item{SNPs}{number of SNPs within the window.}
}
\references{
Waples R. S.: A generalized approach for estimating effective population size from temporal changes in allele frequency, \emph{Genetics} \bold{1989}, 121, 379–391.

Jorde P. E. and Ryman N.: Unbiased estimator for genetic drift and effective population size, \emph{Genetics} \bold{2007}, 177 927–935.

Jónás A., Taus T., Kosiol C., Schlötterer C. & Futschik A.: Estimating effective population size from temporal allele frequency changes in experimental evolution, manuscript in preparation.
}
\author{
Thomas Taus
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{estimateNe}} and \code{\link{read.sync}}
}
\examples{
#MISSING: use example dataset (returned after reading small sync file) and plot some results (barplot & genome-wide Ne)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
