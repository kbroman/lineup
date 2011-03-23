\name{disteg}
\alias{disteg}

\title{Calculate distance between two gene expression data sets}

\description{
  Calculate a distance between all pairs of individuals for two gene
  expression data sets
}

\usage{
disteg(cross, pheno, pmark, min.genoprob=0.99,
       k=20, min.classprob=0.8, repeatKNN=TRUE,
       max.selfd=0.3, phenolabel="phenotype",
       weights, weightByLinkage=TRUE,
       map.function=c("haldane", "kosambi", "c-f", "morgan"),
       standardize=TRUE, verbose=TRUE)

}

\arguments{
 \item{cross}{An object of class \code{"cross"} containing data for a
 QTL experiment.  See the help file for \code{\link[qtl]{read.cross}} in
 the R/qtl package (\url{http://www.rqtl.org}).  There must be a
 phenotype named \code{"id"} or \code{"ID"} that contains the
 individual identifiers.}
 \item{pheno}{A data frame of phenotypes (generally gene expression
 data), stored as individuals x phenotypes.  The row names must contain
 individual identifiers.}
 \item{pmark}{Pseudomarkers that are closest to the genes in
 \code{pheno}, as output by \code{\link{find.gene.pseudomarker}}.}
 \item{min.genoprob}{Threshold on genotype probabilities; if maximum
 probability is less than this, observed genotype taken as \code{NA}.}
 \item{k}{Number of nearest neighbors to consider in forming a k-nearest
 neighbor classifier.}
 \item{min.classprob}{Minimum proportion of neighbors with a common
 class to make a class prediction.}
 \item{repeatKNN}{If TRUE, repeat k-nearest neighbor a second time,
 after omitting individuals who seem to not be self-self matches}
 \item{max.selfd}{Min distance from self (as proportion of mismatches
 between observed and predicted eQTL genotypes) to be excluded from the
 second round of k-nearest neighbor.} 
 \item{phenolabel}{Label for expression phenotypes to place in the
 output distance matrix.}
 \item{weights}{If not missing, use these as weights on the genes in the
calculation of inter-individual distances (actually, since we're working
with eQTL rather than transcripts, we take the average weight for the
multiple transcripts that correspond to a give eQTL},
 \item{weightByLinkage}{If TRUE, weight the eQTL to account for their
relative positions (for example, two tightly linked eQTL would each
count about 1/2 of an isolated eQTL)},
 \item{map.function}{Used if \code{weightByLinkage} is TRUE},
 \item{standardize}{If TRUE, standardize the gene columns so that they
 have mean 0 and SD 1.}
 \item{verbose}{if TRUE, give verbose output.}
}

\details{
We consider the expression phenotypes in batches, by which pseudomarker
they are closest to.  For each batch, we pull the genotype probabilities
at the corresponding pseudomarker and use the individuals that are in
common between \code{cross} and \code{pheno} and whose maximum
genotype probability is above \code{min.genoprob}, to form a classifier
of eQTL genotype from expression values, using k-nearest neighbor (the
function \code{\link[class]{knn}}). The classifier is applied to all
individuals with expression data, to give a predicted eQTL genotype.
(If the proportion of the k nearest neighbors with a common class is
less than \code{min.classprob}, the predicted eQTL genotype is left as
\code{NA}.)

If \code{repeatKNN} is TRUE, we repeat the construction of the k-nearest
neighbor classifier after first omitting individuals whose proportion of
mismatches between observed and inferred eQTL genotypes is greater than
\code{max.selfd}.

Finally, we calculate the distance between the observed eQTL genotypes
for each individual in \code{cross} and the inferred eQTL genotypes for
each individual in \code{pheno}, as the proportion of mismatches between
the observed and inferred eQTL genotypes.

[Explain the use of the weights.]
}

\value{
A matrix with \code{nind(cross)} rows and \code{nrow(pheno)} columns,
containing the distances.  The individual IDs are in the row and column
names.  The matrix is assigned class \code{"lineupdist"}.

The names of the genes that were used to construct the classifier are
saved in an attribute \code{"retained"}.  

The observed and inferred eQTL genotypes are saved as attributes
\code{"obsg"} and \code{"infg"}.  

The denominators of the proportions that form the inter-individual
distances are in the attribute \code{"denom"}.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

%\examples{
%}

\seealso{ \code{\link{distee}}, \code{\link{summary.lineupdist}},
  \code{\link{pulldiag}}, \code{\link{omitdiag}},
  \code{\link{findCommonInd}}, \code{\link{find.gene.pseudomarker}},
  \code{\link{calc.locallod}}, \code{\link{plot.lineupdist}} }

\keyword{utilities}
