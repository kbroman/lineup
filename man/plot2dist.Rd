\name{plot2dist}
\alias{plot2dist}

\title{Plot two sets of inter-individual distances against one another}

\description{
  Plot two sets of inter-individual distances against one another,
  colored by self and non-self distances.
}

\usage{
plot2dist(d1, d2, hirow, hicol, xlab, ylab,
         smoothScatter=TRUE, \dots)
}

\arguments{
 \item{d1}{Output of \code{\link{distee}}.}
 \item{d2}{Output of \code{\link{distee}}.}
 \item{hirow}{Names of rows to highlight in green.}
 \item{hicol}{Names of columns to highlight in orange.}
 \item{xlab}{X-axis label (optional)}
 \item{ylab}{Y-axis label (optional)}
 \item{smoothScatter}{If TRUE, plot non-self distances with
 \code{\link[graphics]{smoothScatter}}; if FALSE, use
 \code{\link[graphics]{plot}}.} 
 \item{\dots}{Passed to \code{\link[graphics]{plot}} and
 \code{\link[graphics]{points}}.} 
}

%\details{
%}

\value{
None.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

%\examples{
%}

\seealso{ \code{\link{pulldiag}}, \code{\link{distee}}, \code{\link{summary.lineupdist}} }

\keyword{graphics}
