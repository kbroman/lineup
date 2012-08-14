\name{disteg}
\alias{disteg}

\title{Calculate distance between two gene expression data sets}

\description{
  Calculate a distance between all pairs of individuals for two gene
  expression data sets
}

\usage{
disteg(cross, pheno, pmark, min.genoprob=0.99,
       k=20, min.classprob=0.8, classprob2drop=1, repeatKNN=TRUE,
       max.selfd=0.3, phenolabel="phenotype",
       weightByLinkage=FALSE,
       map.function=c("haldane", "kosambi", "c-f", "morgan"),
       verbose=TRUE)
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
 \item{classprob2drop}{If an individual is inferred to have a genotype
 mismatch with classprob > this value, treat as an outlier and drop from
 the analysis and then repeat the KNN construction without it.}
 \item{repeatKNN}{If TRUE, repeat k-nearest neighbor a second time,
 after omitting individuals who seem to not be self-self matches}
 \item{max.selfd}{Min distance from self (as proportion of mismatches
 between observed and predicted eQTL genotypes) to be excluded from the
 second round of k-nearest neighbor.} 
 \item{phenolabel}{Label for expression phenotypes to place in the
 output distance matrix.}
 \item{weightByLinkage}{If TRUE, weight the eQTL to account for their
relative positions (for example, two tightly linked eQTL would each
count about 1/2 of an isolated eQTL)}
 \item{map.function}{Used if \code{weightByLinkage} is TRUE}
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

If \code{weightByLinkage} is \code{TRUE}, we use weights on the mismatch
proportions for the various eQTL, taking into account their linkage.
Two tightly linked eQTL will each be given half the weight of a single
isolated eQTL.
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

\examples{
##############################
# simulate an eQTL data set
##############################
# genetic map
L <- seq(120, length=8, by=-10)
map <- sim.map(L, n.mar=L/10+1, include.x=FALSE, eq.spacing=TRUE)

# physical map: make all intervals 2x longer
pmap <- rescalemap(map, 2)

# arbitrary locations of 40 local eQTL
thepos <- unlist(map)
theppos <- unlist(pmap)
thechr <- rep(seq(along=map), sapply(map, length))
eqtl.loc <- sort(sample(seq(along=thepos), 40))

x <- sim.cross(map, n.ind=250, type="f2",
               model=cbind(thechr[eqtl.loc], thepos[eqtl.loc], 0, 0))
x$pheno$id <- factor(paste("Mouse", 1:250, sep=""))

# first 20 have eQTL with huge effects
# second 20 have essentially no effect
edata <- cbind((x$qtlgeno[,1:20] - 2)*10+rnorm(prod(dim(x$qtlgeno[,1:20]))),
               (x$qtlgeno[,21:40] - 2)*0.1+rnorm(prod(dim(x$qtlgeno[,21:40]))))
dimnames(edata) <- list(x$pheno$id, paste("e", 1:ncol(edata), sep=""))

# gene locations
theloc <- data.frame(chr=thechr[eqtl.loc], pos=theppos[eqtl.loc])
rownames(theloc) <- colnames(edata)

# mix up 5 individuals in expression data
edata[1:3,] <- edata[c(2,3,1),]
edata[4:5,] <- edata[5:4,]

##############################
# now, the start of the analysis
##############################
x <- calc.genoprob(x, step=1)

# find nearest pseudomarkers
pmark <- find.gene.pseudomarker(x, pmap, theloc, "prob")

# calculate LOD score for local eQTL
locallod <- calc.locallod(x, edata, pmark)

# take those with LOD > 100 [which will be the first 20]
edatasub <- edata[,locallod>100,drop=FALSE]

# calculate distance between individuals
#     (prop'n mismatches between obs and inferred eQTL geno)
d <- disteg(x, edatasub, pmark)

# plot distances
plot(d)

# summary of apparent mix-ups
summary(d)

# plot of classifier for first eQTL
plotEGclass(d)
}

\seealso{ \code{\link{distee}}, \code{\link{summary.lineupdist}},
  \code{\link{pulldiag}}, \code{\link{omitdiag}},
  \code{\link{findCommonID}}, \code{\link{find.gene.pseudomarker}},
  \code{\link{calc.locallod}}, \code{\link{plot.lineupdist}},
  \code{\link[class]{knn}}, \code{\link{plotEGclass}} }

\keyword{utilities}
