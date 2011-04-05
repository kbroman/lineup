\name{distee}
\alias{distee}

\title{Calculate distance between two gene expression data sets}

\description{
  Calculate a distance between all pairs of individuals for two gene
  expression data sets
}

\usage{
distee(e1, e2, d.method=c("rmsd", "cor"), labels=c("e1","e2"),
       verbose=TRUE)
}

\arguments{
 \item{e1}{Numeric matrix of gene expression data, as individuals x
 genes.  The row and column names must contain individual and gene
 identifiers.}
 \item{e2}{(Optional) Like \code{e1}.  An appreciable number of
 individuals and genes must be in common.}
 \item{d.method}{Calculate inter-individual distance as RMS difference
 or as correlation.}
 \item{labels}{Two character strings, to use as labels for the two data
 matrices in subsequent output.}
 \item{verbose}{if TRUE, give verbose output.}
}

\details{
We calculate the pairwise distance between all individuals (rows) in
\code{e1} and all individuals in \code{e2}.  This distance is either the
RMS difference (\code{d.method="rmsd"}) or the correlation
(\code{d.method="cor"}).
}

\value{
A matrix with \code{nrow(e1)} rows and \code{nrow(e2)} columns,
containing the distances.  The individual IDs are in the row and column
names.  The matrix is assigned class \code{"lineupdist"}.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

\examples{
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

\seealso{ \code{\link{pulldiag}}, \code{\link{omitdiag}},
\code{\link{summary.lineupdist}}, \code{\link{plot2dist}},
\code{\link{disteg}}, \code{\link{corbetw2mat}} }

\keyword{utilities}
