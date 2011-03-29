\name{find.gene.pseudomarker}
\alias{find.gene.pseudomarker}

\title{Find nearest peudomarker to each gene}

\description{
  Pull out the pseudomarker that is closest to the position of each of a
  series of genes.
}

\usage{
find.gene.pseudomarker(cross, pmap, geneloc, where=c("prob", "draws"))
}

\arguments{
 \item{cross}{An object of class \code{"cross"} containing data for a
 QTL experiment.  See the help file for \code{\link[qtl]{read.cross}} in
 the R/qtl package (\url{http://www.rqtl.org}).}
 \item{pmap}{A physical map of the markers in \code{cross}, with
 locations in Mbp.  This is a list whose components are the marker
 locations on each chromosome.}
 \item{geneloc}{A data frame specifying the physical locations of the
 genes.  There should be two columns, \code{chr} for chromosome and
 \code{pos} for position in Mbp.  The rownames should indicate the gene
 names.}
 \item{where}{Indicates whether to pull pseudomarkers from the genotype
 probabilities (produced by \code{\link[qtl]{calc.genoprob}}) or from
 the imputed genotypes (produced by \code{\link[qtl]{sim.geno}}).}
}

\details{
  We first convert positions (by interpolation) from those contained
  within \code{cross} to physical coordinates contained in \code{pmap}.
  We then use \code{\link[qtl]{find.pseudomarker}} to identify the
  closest pseudomarker to each gene location.

  We also include the positions of the pseudomarkers, and we print a
  warning message if pseudomarkers are > 2 Mbp from the respective gene.
}

\value{
A data frame with columns \code{chr} (the chromosome) and \code{pmark}
(the name of the pseudomarker).  The third column \code{pos} contains
the Mbp position of the pseudomarker.  The final column is the signed
distance between the gene and the pseudomarker.  The rownames indicate
the gene names.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

\examples{
\dontrun{
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
}

\seealso{ \code{\link[qtl]{find.pseudomarker}},
  \code{\link[qtl]{find.pseudomarkerpos}}, \code{\link{plotEGclass}},
  \code{\link{disteg}}, \code{\link{calc.locallod}} }

\keyword{utilities}
