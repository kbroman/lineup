\name{distee}
\alias{distee}

\title{Calculate distance between two gene expression data sets}

\description{
  Calculate a distance between all pairs of individuals for two gene
  expression data sets
}

\usage{
distee(e1, e2, cor.threshold, n.col, d.method=c("rmsd", "cor"),
       transpose=FALSE, labels=c("e1","e2"), verbose=TRUE)
}

\arguments{
 \item{e1}{Numeric matrix of gene expression data, as individuals x
 genes (or, if \code{transpose=TRUE}, genes x individuals).  The row and
 column names must contain individual and gene identifiers.}
 \item{e2}{Like \code{e1}.  An appreciable number of individuals and
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
}

\value{
A matrix with \code{nrow(e1)} rows and \code{nrow(e2)} columns,
containing the distances.  The individual IDs are in the row and column
names.  The matrix is assigned class \code{"lineupdist"}.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

%\examples{
%}

\seealso{ \code{\link{pulldiag}}, \code{\link{summary.lineupdist}} }

\keyword{utilities}
