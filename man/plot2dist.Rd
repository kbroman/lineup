\name{plot2dist}
\alias{plot2dist}

\title{Plot two sets of inter-individual distances against one another}

\description{
  Plot two sets of inter-individual distances against one another,
  colored by self and non-self distances.
}

\usage{
plot2dist(d1, d2, hirow, hicol, xlab, ylab, smoothScatter=FALSE,
         colself="black", colnonself="gray", colhirow="green", colhicol="orange", \dots)
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
 \item{colself}{Color to use for the self-self points.  If NULL,
 these aren't plotted.}
 \item{colnonself}{Color to use for the non-self points.  If NULL,
 these aren't plotted.}
 \item{colhirow}{Color to use for the \code{hirow} points.  If NULL,
 these aren't plotted.}
 \item{colhicol}{Color to use for the \code{hicol} points.  If NULL,
 these aren't plotted.}
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

\seealso{ \code{\link{pulldiag}}, \code{\link{distee}}, \code{\link{summary.lineupdist}} }

\keyword{graphics}
