\name{scaleAF}
\alias{scaleAF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scale Allele Frequencies
}
\description{
Allele frequencies are scaled by logit- or arcsin-transformation.
}
\usage{
scaleAF(af, method = c("logit", "asin"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{af}{
numeric vector of allele frequencies [0,1] that should be transformed
}
  \item{method}{
character string determining the scaling method. Possible values are \code{"logit"} (default) and \code{"asin"}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Numeric vector of transformed allele frequencies
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
