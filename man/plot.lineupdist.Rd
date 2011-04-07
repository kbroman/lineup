\name{plot.lineupdist}
\alias{plot.lineupdist}

\title{Plot summary of inter-individual distances}

\description{
  Plot histograms of self-self and self-nonself distances from a
  distance matrix calculated by \code{\link{distee}} or \code{\link{disteg}}.
}

\usage{
\method{plot}{lineupdist}(x, breaks, add.rug=TRUE, what=c("both", "ss", "sn"), \dots)
}

\arguments{
 \item{x}{Output of \code{\link{distee}} or \code{\link{disteg}}.}
 \item{breaks}{Optional vector of breaks, passed to
   \code{\link[graphics]{hist}}, though if it is length 1, we interpret
 it as the number of breaks and ensure that both histograms use the same
 set of breaks.}
 \item{add.rug}{If true, also include \code{\link[graphics]{rug}} below histograms.}
 \item{what}{Indicates whether to plot both self-self and self-nonself
 distances (or correlations) or just one or the other.  (\code{"ss"}
 indicates self-self and \code{"sn"} indicates self-nonself.)}
 \item{\dots}{Ignored at this point.}
}

\details{
  We call \code{\link{pulldiag}} and \code{\link{omitdiag}} to get the
  self-self and self-nonself distances.

  If all of the self-self distances are missing, we plot just the
  self-nonself distances.
}

\value{
None.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

\examples{
\dontrun{
# simulate MVN, 100 individuals, 40 measurements (of which 20 are just noise)
V <- matrix(0.3, ncol=20, nrow=20) + diag(rep(0.5, 20)) 
D <- chol(V)
z <- matrix(rnorm(20*100), ncol=20) %*% D

# create two data matrices as z + noise
x <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))
y <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))

# permute some rows
x[51:53,] <- x[c(52,53,51),]
y[41:42,] <- y[42:41,]

# add column and row names
dimnames(x) <- dimnames(y) <- list(paste("ind", 1:100, sep=""),
                                   paste("gene", 1:40, sep=""))

# calculate correlations between cols of x and cols of y
thecor <- corbetw2mat(x, y)

# subset x and y, taking only columns with corr > 0.75
xs <- x[,thecor > 0.8]
ys <- y[,thecor > 0.8]

# calculate distance (using "RMS difference" as a measure)
d1 <- distee(xs, ys, d.method="rmsd", labels=c("x","y"))

# calculate distance (using "correlation" as a measure...really similarity)
d2 <- distee(xs, ys, d.method="cor", labels=c("x", "y"))

# pull out the smallest 8 self-self correlations
sort(pulldiag(d2))[1:8]

# summary of results
summary(d1)
summary(d2)

# order to put matches together
summary(d2, reorder="alignmatches")

# plot histograms of RMS distances
plot(d1)

# plot histograms of correlations
plot(d2)

# plot distances against one another
plot2dist(d1, d2)
}
}

\seealso{ \code{\link{pulldiag}}, \code{\link{distee}}, \code{\link{plot2dist}} }

\keyword{utilities}
