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

\examples{
\dontrun{
# simulate MVN, 100 individuals, 20 measurements
V <- matrix(0.3, ncol=20, nrow=20) + diag(rep(0.5, 20)) 
D <- chol(V)
z <- matrix(rnorm(20*100), ncol=20) %*% D

# create two data matrices as z + noise
x <- z + rnorm(20*100, 0, 0.2)
y <- z + rnorm(20*100, 0, 0.2)

# permute some rows
x[51:53,] <- x[c(52,53,51),]
y[41:42,] <- y[42:41,]

# add column and row names
dimnames(x) <- dimnames(y) <- list(paste("ind", 1:100, sep=""),
                                   paste("gene", 1:20, sep=""))

# calculate distance (using "RMS difference" as a measure)
d1 <- distee(x, y, cor.threshold=0.8)

# calculate distance (using "correlation" as a measure...really similarity)
d2 <- distee(x, y, cor.threshold=0.8, d.method="cor")

# pull out the smallest 8 self-self correlations
sort(pulldiag(d2))[1:8]

# summary of results
summary(d1)
summary(d2)

# plot distances against one another
plot2dist(d1, d2)
}
}

\seealso{ \code{\link{pulldiag}}, \code{\link{distee}}, \code{\link{plot2dist}} }

\keyword{utilities}
