\name{sync}
\alias{sync-class}
\alias{is.sync}
\alias{merge,sync,sync-method}



\title{Sync Object - Creation and Usage}

\description{\code{sync} is an S4 class that stores all information of a synchronized file generated e.g. by the \href{http://sourceforge.net/projects/popoolation2/}{Popoolation2} software. }

\usage{
is.sync(x)

\S4method{merge}{sync,sync}(x, y, ...)
}

\arguments{
  \item{x, y}{
any \R object for \code{is.sync}; an object inheriting from class "sync" otherwise.
}
  \item{...}{
not used at the moment.
  }
}

\details{
\code{is.sync} checks if an object is a \code{sync} object.
}

\value{
\code{is.sync} returns \code{TRUE} if its argument is a \code{sync} object and \code{FALSE} otherwise.
}

\author{
Thomas Taus
}

\seealso{
\code{\link{read.sync}}, \code{\link{alleles}}, \code{\link{polarization}}, \code{\link{af.traj}}, \code{\link{af}} and \code{\link{coverage}}
}
