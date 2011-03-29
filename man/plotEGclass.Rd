\name{plotEGclass}
\alias{plotEGclass}

\title{Plot classifier of eQTL genotype from expression data}

\description{
  Diagnostic plot of one of the eQTL classifiers from the results of
  \code{\link{disteg}}: generally expression phenotype against observed
  eQTL genotype, colored by inferred eQTL genotype.
}

\usage{
plotEGclass(d, eqtl=1, outercol="inferred", innercol="observed", 
            thecolors=c("blue","green","red","orange"), \dots)
}

\arguments{
  \item{d}{Output of \code{\link{disteg}}.}
  \item{eqtl}{Numeric index or a character vector (of the form
    "1@102.35") indicating the eQTL to consider.}
  \item{outercol}{Indicates how to color the outer edge of the points:
  \code{"observed"} indicates to color based on observed genotypes;
  \code{"inferred"} indicates to color based on inferred genotypes;
  otherwise, give a color.}
  \item{innercol}{Like \code{outercol}, but indicating the interior of
    the points.}
  \item{thecolors}{The colors to use in the plot.  The last element
  (after the number of genotypes) indicates the color to use for missing
  values.}
 \item{\dots}{Passed to \code{\link[graphics]{plot}} and
 \code{\link[graphics]{points}}.} 
}

\details{
  The function produces a diagnostic plot for studying one of the
  k-nearest neighbor classifiers underlying the output from
  \code{\link{disteg}}. 

  In the case of one expression phenotype attached to the selected eQTL, 
  the plot is a dot plot of gene expression against observed eQTL genotype.

  In the case of two expression phenotypes, the plot is a scatterplot of
  the two expression phenotypes against each other.

  In the case of more than two expression phenotypes, we use
  \code{\link[graphics]{pairs}} to produce a matrix of scatterplots.
}

\value{
None.
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

\seealso{ \code{\link{disteg}}, \code{\link{plot.lineupdist}},
  \code{\link{plot2dist}}, \code{\link[class]{knn}} }

\keyword{graphics}
