\name{pulldiag}
\alias{pulldiag}

\title{Pull out the diagonal from a distance matrix}

\description{
  Pull out the diagonal from a distance matrix calculated by
  \code{\link{distee}} (that is, self-self distances).
}

\usage{
pulldiag(d)
}

\arguments{
 \item{d}{A distance matrix calculated by \code{\link{distee}}.}
}

\details{
We use the row and column names to identify which entries are self-self.
}

\value{
A vector with the self-self distances.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

%\examples{
%}

\seealso{ \code{\link{distee}} }

\keyword{utilities}
