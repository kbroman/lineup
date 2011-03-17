\name{summary.lineupdist}
\alias{summary.lineupdist}

\title{Summarize inter-individual distances}

\description{
  Summarize the results of \code{\link{distee}}, with inter-individual
  distances between two sets of gene expression data.
}

\usage{
\method{summary}{lineupdist}(object, dropmatches=TRUE, reorder=TRUE,
                  \dots)
}

\arguments{
 \item{object}{Output of \code{\link{distee}}.}
 \item{dropmatches}{If TRUE, omit rows for which an individual's best
 match is itself.}
 \item{reorder}{If TRUE, reorder rows by increasing distance (or
 decreasing correlation) to the best match and then by decreasing
 distance (or decreasing correlation) to self.}
 \item{\dots}{Passed to \code{\link[base]{print.data.frame}}.}
}

%\details{
%}

\value{
A list with two components: the distances summarized by row and the
distances summarized by column.

For each individual, we calculate the minimum distance to others,
next-smallest distance, the self-self distance, the mean and SD of the
distances to others, and finally indicate the individual (or
individuals) that is closest.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

%\examples{
%}

\seealso{ \code{\link{pulldiag}}, \code{\link{distee}}, \code{\link{plot2dist}} }

\keyword{utilities}
