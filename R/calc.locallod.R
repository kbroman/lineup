## calc.locallod.R
## Karl W Broman

# calc.locallod
#
#' Calculate LOD score at physical position of each gene
#'
#' For gene expression data with physical positions of the genes, calculate the
#' LOD score at those positions to assess evidence for local eQTL.
#'
#' \code{cross} and \code{pheno} must contain exactly the same individuals in
#' the same order.  (Use \code{\link{findCommonID}} to line them up.)
#'
#' We consider the expression phenotypes in batches: those whose closest
#' pseudomarker is the same.
#'
#' We use Haley-Knott regression to calculate the LOD scores.
#'
#' Actually, we use a bit of a contortion of the data to force the
#' \code{\link[qtl]{scanone}} function in R/qtl to calculate the LOD score at a
#' single position.
#'
#' We omit any transcripts that map to the X chromosome; we can only handle
#' autosomal loci for now.
#'
#' @param cross An object of class \code{"cross"} containing data for a QTL
#' experiment.  See the help file for \code{\link[qtl]{read.cross}} in the
#' R/qtl package (\url{http://www.rqtl.org}).  There must be a phenotype named
#' \code{"id"} or \code{"ID"} that contains the individual identifiers.
#' @param pheno A data frame of phenotypes (generally gene expression data),
#' stored as individuals x phenotypes.  The row names must contain individual
#' identifiers.
#' @param pmark Pseudomarkers that are closest to the genes in \code{pheno}, as
#' output by \code{\link{find.gene.pseudomarker}}.
#' @param addcovar Additive covariates passed to \code{\link{scanone}}.
#' @param intcovar Interactive covariates passed to \code{\link{scanone}}.
#' @param verbose If TRUE, print tracing information.
#' @return A vector of LOD scores.  The names indicate the gene names (rows in
#' \code{pheno}).
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{find.gene.pseudomarker}}, \code{\link{plotEGclass}},
#' \code{\link{findCommonID}}, \code{\link{disteg}}
#' @keywords utilities
#' @examples
#'
#' \dontrun{
#' ##############################
#' # simulate an eQTL data set
#' ##############################
#' # genetic map
#' L <- seq(120, length=8, by=-10)
#' map <- sim.map(L, n.mar=L/10+1, include.x=FALSE, eq.spacing=TRUE)
#'
#' # physical map: make all intervals 2x longer
#' pmap <- rescalemap(map, 2)
#'
#' # arbitrary locations of 40 local eQTL
#' thepos <- unlist(map)
#' theppos <- unlist(pmap)
#' thechr <- rep(seq(along=map), sapply(map, length))
#' eqtl.loc <- sort(sample(seq(along=thepos), 40))
#'
#' x <- sim.cross(map, n.ind=250, type="f2",
#'                model=cbind(thechr[eqtl.loc], thepos[eqtl.loc], 0, 0))
#' x$pheno$id <- factor(paste("Mouse", 1:250, sep=""))
#'
#' # first 20 have eQTL with huge effects
#' # second 20 have essentially no effect
#' edata <- cbind((x$qtlgeno[,1:20] - 2)*10+rnorm(prod(dim(x$qtlgeno[,1:20]))),
#'                (x$qtlgeno[,21:40] - 2)*0.1+rnorm(prod(dim(x$qtlgeno[,21:40]))))
#' dimnames(edata) <- list(x$pheno$id, paste("e", 1:ncol(edata), sep=""))
#'
#' # gene locations
#' theloc <- data.frame(chr=thechr[eqtl.loc], pos=theppos[eqtl.loc])
#' rownames(theloc) <- colnames(edata)
#'
#' # mix up 5 individuals in expression data
#' edata[1:3,] <- edata[c(2,3,1),]
#' edata[4:5,] <- edata[5:4,]
#'
#' ##############################
#' # now, the start of the analysis
#' ##############################
#' x <- calc.genoprob(x, step=1)
#'
#' # find nearest pseudomarkers
#' pmark <- find.gene.pseudomarker(x, pmap, theloc, "prob")
#'
#' # calculate LOD score for local eQTL
#' locallod <- calc.locallod(x, edata, pmark)
#'
#' # take those with LOD > 100 [which will be the first 20]
#' edatasub <- edata[,locallod>100,drop=FALSE]
#'
#' # calculate distance between individuals
#' #     (prop'n mismatches between obs and inferred eQTL geno)
#' d <- disteg(x, edatasub, pmark)
#'
#' # plot distances
#' plot(d)
#'
#' # summary of apparent mix-ups
#' summary(d)
#'
#' # plot of classifier for first eQTL
#' plotEGclass(d)
#' }
#'
#' @importFrom qtl nind scanone
#' @export calc.locallod
calc.locallod <-
    function(cross, pheno, pmark, addcovar=NULL, intcovar=NULL, verbose=TRUE)
{
    if(any(pmark$chr == "X")) {
        warning("Dropping X chr loci; we can only handle autosomes for now.")
        pmark <- pmark[pmark$chr != "X",]
    }

    if(nind(cross) != nrow(pheno))
        stop("cross and pheno have incompatible numbers of individuals.")
    m <- match(colnames(pheno), rownames(pmark))
    if(any(is.na(m))) pheno <- pheno[,!is.na(m),drop=FALSE]
    m <- match(rownames(pmark), colnames(pheno))
    if(any(is.na(m))) pmark <- pmark[!is.na(m),,drop=FALSE]
    pheno <- pheno[,match(rownames(pmark), colnames(pheno)),drop=FALSE]

    cpmark <- paste(pmark$chr, pmark$pmark, sep=":")
    upmark <- unique(cpmark)

    lod <- rep(NA, ncol(pheno))
    temp <- cross
    n.ind <- nind(cross)

    # loop over unique pseudomarkers
    for(i in seq(along=upmark)) {
        if(verbose && i==round(i,-2)) cat(i, "of", length(upmark), "\n")
        wh <- which(cpmark == upmark[i])

        y <- pheno[,wh,drop=FALSE]
        gp <- cross$geno[[pmark$chr[wh[1]]]]$prob[,pmark$pmark[wh[1]],,drop=FALSE]

        # create dummy cross
        temp$pheno <- cbind(y, cross$pheno)
        temp$geno <- list("1"=list(data=cbind("m1"=rep(1, n.ind)), map=c("m1"=1),prob=gp))
        attr(temp$geno[[1]]$prob, "map") <- c("m1"=1)
        class(temp$geno[[1]]) <- "A"
        lod[wh] <- unlist(scanone(temp, method="hk", addcovar=addcovar, intcovar=intcovar,
                                  pheno.col=1:ncol(y))[-(1:2)])
    }

    names(lod) <- colnames(pheno)
    lod
}
