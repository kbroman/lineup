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
the Mbp position of the pseudomarker.  The rownames indicate the gene
names.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

%\examples{
%}

\seealso{ \code{\link[qtl]{find.pseudomarker}},
  \code{\link[qtl]{find.pseudomarkerpos}}, 
  \code{\link{disteg}} \code{\link{calc.locallod}} }

\keyword{utilities}
