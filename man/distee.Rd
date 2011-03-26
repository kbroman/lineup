\name{distee}
\alias{distee}

\title{Calculate distance between two gene expression data sets}

\description{
  Calculate a distance between all pairs of individuals for two gene
  expression data sets
}

\usage{
distee(e1, e2, cor.threshold, n.col, d.method=c("rmsd", "cor"),
       transpose=FALSE, labels=c("e1","e2"), weightByCorr=FALSE,
       standardize=FALSE, verbose=TRUE)
}

\arguments{
 \item{e1}{Numeric matrix of gene expression data, as individuals x
 genes (or, if \code{transpose=TRUE}, genes x individuals).  The row and
 column names must contain individual and gene identifiers.}
 \item{e2}{(Optional) Like \code{e1}.  An appreciable number of individuals and
 genes must be in common.}
 \item{cor.threshold}{threshold on between-gene correlation; use genes
 with between-condition correlation above this to calculate the
 distance.}
 \item{n.col}{if \code{cor.threshold} not given, instead provide the
 number of genes to use (specified by this argument).  This number will
 be chosen, in decreasing order of correlation.}
 \item{d.method}{Calculate inter-individual distance as RMS difference
 or as correlation.}
 \item{transpose}{if TRUE, \code{e1} and \code{e2} are input as genes x individuals
 and need to be transposed.}
 \item{labels}{Two character strings, to use as labels for the two data
 matrices in subsequent output.}
 \item{weightByCorr}{If TRUE, in the calculate of inter-individual
 distances, weight genes by their cross-tissue correlation.}
 \item{standardize}{If TRUE, standardize the gene columns so that they
 have mean 0 and SD 1.}
 \item{verbose}{if TRUE, give verbose output.}
}

\details{
We first line up the data matrices, dropping genes that are not in
common between the two.  Each column (a gene) is standardized to have
mean 0 and SD 1.

Using individuals in common between the two data sets, we calculate, for
each gene, the correlation between the two conditions across
individuals.  If \code{cor.threshold} is given, genes with correlation
\eqn{\ge}{>=} \code{cor.threshold} are retained.  If \code{n.col} is
given, the top \code{n.col} genes, in decreasing order of correlation,
are retained.

Finally, we calculate the pairwise distance between all individuals
(rows) in \code{e1} and all individuals in \code{e2}.  This distance is
either the RMS difference (\code{d.method="rmsd"}) or the correlation
(\code{d.method="cor"}).  

If \code{e2} is not included, within-tissue comparisons will be made.
In this case, genes will have to be pre-selected (e.g., as those
selected in a comparison of one tissue to another).

[Explain the use of the weights.]
}

\value{
A matrix with \code{nrow(e1)} rows and \code{nrow(e2)} columns,
containing the distances.  The individual IDs are in the row and column
names.  The matrix is assigned class \code{"lineupdist"}.

The names of the genes that are retained are saved in an attribute
\code{"retained"}.
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

# calculate distance (using "RMS difference" as a measure)
d1 <- distee(x, y, cor.threshold=0.8)

# calculate distance (using "correlation" as a measure...really similarity)
d2 <- distee(x, y, cor.threshold=0.8, d.method="cor")

# pull out the smallest 8 self-self correlations
sort(pulldiag(d2))[1:8]

# summary of results
summary(d1)
summary(d2)

# plot histograms of RMS distances
plot(d1)

# plot histograms of correlations
plot(d2)

# plot distances against one another
plot2dist(d1, d2)
}

\seealso{ \code{\link{pulldiag}}, \code{\link{omitdiag}},
\code{\link{summary.lineupdist}}, \code{\link{plot2dist}},
\code{\link{disteg}}, \code{\link{corbetw2mat}}, \code{\link{fscale}} }

\keyword{utilities}
