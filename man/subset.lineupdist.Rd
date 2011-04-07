\name{subset.lineupdist}
\alias{subset.lineupdist}
\alias{[.lineupdist}

\title{Subsetting distance matrix}

\description{
  Pull out a specified set of rows and columns from a distance matrix
  calculated by \code{\link{distee}} or \code{\link{disteg}}.
}

\usage{
\method{subset}{lineupdist}(x, rows, cols, \dots)
\method{[}{lineupdist}(x, rows, cols)
}

\arguments{
 \item{x}{A distance matrix object as obtained from \code{\link{distee}}
 or \code{\link{disteg}}.}
 \item{rows}{Optional vector of selected rows.}
 \item{cols}{Optional vector of selected columns.}
 \item{\dots}{Ignored at this point.}
}

\value{The input distance matrix object, but with only the specified subset
  of the data.}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\examples{
%}

\seealso{ \code{\link{disteg}}, \code{\link{distee}}, \code{\link{pulldiag}} }

\keyword{manip}
