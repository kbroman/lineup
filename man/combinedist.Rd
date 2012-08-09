\name{combinedist}
\alias{combinedist}

\title{Combine distance matrices into a single such}

\description{
  Combine multiple distance matrices into a single distance matrix
  providing an overall summary
}

\usage{
combinedist(\dots, method=c("median", "mean"))
}

\arguments{
 \item{\dots}{Set of distance matrices, as calculated by
 \code{\link{distee}} or \code{\link{disteg}}.}
 \item{method}{Indicates whether to summarize using the median or the mean.}
}

\details{
  The row and column names of the input distance matrices define the
  individual IDs.

  If the input distance matrices all have an attribute \code{"denom"}
  (for denominator) and \code{method="mean"}, we use a weighted mean,
  weighted by the denominators.  This could be used to calculate an
  overall proportion.
}

\value{
A distance matrix, with class \code{"lineupdist"}. 
The individual IDs are in the row and column
names.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

\examples{
# simulate MVN, 100 individuals, 40 measurements (of which 20 are just noise)
V <- matrix(0.3, ncol=20, nrow=20) + diag(rep(0.5, 20)) 
D <- chol(V)
z <- matrix(rnorm(20*100), ncol=20) %*% D

# create three data matrices as z + noise
x <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))
y <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))
w <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))

# permute some rows of x
x[51:53,] <- x[c(52,53,51),]

# add column and row names
dimnames(x) <- dimnames(y) <- dimnames(w) <- 
   list(paste("ind", 1:100, sep=""), paste("gene", 1:40, sep=""))

# calculate correlations between cols of x and of the other two matrices
corxy <- corbetw2mat(x, y)
corxw <- corbetw2mat(x, w)

# using columns with corr > 0.75, 
# calculate distance (using "correlation" as a measure...really similarity)
dxy <- distee(x[,corxy>0.75], y[,corxy>0.75], d.method="cor", labels=c("x", "y"))
dxw <- distee(x[,corxw>0.75], w[,corxw>0.75], d.method="cor", labels=c("x", "w"))

d <- combinedist(dxy, dxw)

summary(d)
}

\seealso{ \code{\link{distee}}, \code{\link{disteg}},
\code{\link{summary.lineupdist}} }

\keyword{utilities}
