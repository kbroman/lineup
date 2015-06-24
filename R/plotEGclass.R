## plotEGclass.R
## Karl W Broman

# plotEGclass
#
#' Plot classifier of eQTL genotype from expression data
#'
#' Diagnostic plot of one of the eQTL classifiers from the results of
#' \code{\link{disteg}}: generally expression phenotype against observed eQTL
#' genotype, colored by inferred eQTL genotype.
#'
#' The function produces a diagnostic plot for studying one of the k-nearest
#' neighbor classifiers underlying the output from \code{\link{disteg}}.
#'
#' In the case of one expression phenotype attached to the selected eQTL, the
#' plot is a dot plot of gene expression against observed eQTL genotype.
#'
#' In the case of two expression phenotypes, the plot is a scatterplot of the
#' two expression phenotypes against each other.
#'
#' In the case of more than two expression phenotypes, we use
#' \code{\link[graphics]{pairs}} to produce a matrix of scatterplots.
#'
#' @param d Output of \code{\link{disteg}}.
#' @param eqtl Numeric index or a character vector (of the form "1@@102.35")
#' indicating the eQTL to consider.
#' @param outercol Indicates how to color the outer edge of the points:
#' \code{"observed"} indicates to color based on observed genotypes;
#' \code{"inferred"} indicates to color based on inferred genotypes; otherwise,
#' give a color.
#' @param innercol Like \code{outercol}, but indicating the interior of the
#' points.
#' @param thecolors The colors to use in the plot.  The last element (after the
#' number of genotypes) indicates the color to use for missing values.
#' @param \dots Passed to \code{\link[graphics]{plot}} and
#' \code{\link[graphics]{points}}.
#' @return None.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{disteg}}, \code{\link{plot.lineupdist}},
#' \code{\link{plot2dist}}, \code{\link[class]{knn}}
#' @keywords graphics
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
#' @importFrom graphics plot par pairs points
#' @export
plotEGclass <-
    function(d, eqtl=1, outercol="inferred", innercol="observed",
             thecolors=c("blue","green","red","orange"), ...)
{
    if(!("eg.lineupdist" %in% class(d)))
        stop("Input d must be as produced by disteg().")

    # inner and outer colors based on...
    choices <- c("observed", "inferred")
    if(!is.na(pm <- pmatch(outercol, choices)))
        outercol <- choices[pm]
    if(!is.na(pm <- pmatch(innercol, choices)))
        innercol <- choices[pm]

    ginf <- attr(d, "infg")
    gobs <- attr(d, "obsg")
    y <- attr(d, "y")
    gnames <- attr(d, "genonames")

    if(length(eqtl) > 1) {
        warning("eqtl should have length 1; using first element.")
        eqtl <- eqtl[1]
    }
    if(is.character(eqtl)) { # name of an eQTL
        eqtlnam <- eqtl
        eqtl <- match(eqtl, colnames(ginf))
        if(is.na(eqtl))
            stop("Can't find eQTL ", eqtl)
    }
    else {
        if(eqtl < 1 || eqtl > ncol(ginf))
            stop("eqtl must be between 1 and ", ncol(ginf))

        eqtlnam <- colnames(ginf)[eqtl]
    }

    # clean up eQTL name
    if(length(grep("@", eqtlnam)) == 1) {
        spl <- unlist(strsplit(eqtlnam, "@"))
        splnum <- unlist(strsplit(as.character(round(as.numeric(spl[2]),2)), "\\."))
        if(length(splnum)==1 || nchar(splnum[2]) == 0)
            splnum[2] <- "00"
        else if(nchar(splnum[2]) == 1)
            splnum[2] <- paste(splnum[2], "0", sep="")
        eqtlnam <- paste(spl[1], " @ ", splnum[1], ".", splnum[2], sep="")
    }

    ids <- findCommonID(rownames(gobs), rownames(ginf))

    ginf <- ginf[,eqtl]
    gobs <- gobs[,eqtl]
    y <- y[[eqtl]]
    hasy <- rowSums(is.na(y)) == 0

    gobs.sub <- gobs[ids$first]
    ginf.sub <- ginf[ids$second]
    y.sub <- y[ids$second,,drop=FALSE]

    ynog <- y[-ids$second,,drop=FALSE]
    ynog <- ynog[rowSums(is.na(ynog))==0,,drop=FALSE]

    if(ncol(y)==1) { # dot plot
        if(nrow(ynog) > 0) gnames <- c(gnames, "NA")

        xjit <- 0.2
        u <- runif(length(gobs.sub), -xjit, xjit)

        # point colors
        obscol <- rep("", length(gobs.sub))
        gobs.sub[is.na(gobs.sub)] <- 4
        for(i in seq(along=gnames))
            obscol[gobs.sub==i] <- thecolors[i]
        infcol <- rep("", length(ginf.sub))
        ginf.sub[is.na(ginf.sub)] <- 4
        for(i in seq(along=gnames))
            infcol[ginf.sub==i] <- thecolors[i]

        col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
        bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

        plot(gobs.sub+u, y.sub, xlab="Observed genotype", ylab=colnames(y),
             main=eqtlnam, xaxt="n", xlim=c(0.5, length(gnames)+0.5),
             pch=21, col=col, bg=bg, ylim=range(y, na.rm=TRUE), type="n", ...)
        axis(side=1, at=seq(along=gnames), gnames)

        # make sure the mismatches are on top
        wh <- !is.na(gobs.sub) & !is.na(ginf.sub) & gobs.sub==ginf.sub
        points((gobs.sub+u)[wh], y.sub[wh],
               pch=21, col=col[wh], bg=bg[wh], ...)
        points((gobs.sub+u)[!wh], y.sub[!wh],
               pch=21, col=col[!wh], bg=bg[!wh], ...)

        if(nrow(ynog) > 0) { # phenotype no genotype
            theinf <- ginf[rownames(ynog)]

            infcol <- rep("", length(nrow(ynog)))
            theinf[is.na(theinf)] <- 4
            for(i in seq(along=gnames))
                infcol[theinf==i] <- thecolors[i]
            obscol <- rep("black", length(infcol))

            col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
            bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

            points(runif(nrow(ynog), 4-xjit, 4+xjit), ynog,
                   pch=21, col=col, bg=bg, ...)
        }

    }
    else if(ncol(y) == 2) { # scatter plot

        if(nrow(ynog) > 0) gnames <- c(gnames, "NA")

        # point colors
        obscol <- rep("", length(gobs.sub))
        gobs.sub[is.na(gobs.sub)] <- 4
        for(i in seq(along=gnames))
            obscol[gobs.sub==i] <- thecolors[i]
        infcol <- rep("", length(ginf.sub))
        ginf.sub[is.na(ginf.sub)] <- 4
        for(i in seq(along=gnames))
            infcol[ginf.sub==i] <- thecolors[i]

        col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
        bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

        plot(y.sub[,1], y.sub[,2], xlab=colnames(y)[1], ylab=colnames(y)[2],
             main=eqtlnam, pch=21, col=col, bg=bg,
             xlim=range(y[,1], na.rm=TRUE), ylim=range(y[,2], na.rm=TRUE), type="n", ...)

        # make sure the mismatches are on top
        wh <- !is.na(gobs.sub) & !is.na(ginf.sub) & gobs.sub==ginf.sub
        points(y.sub[wh,1], y.sub[wh,2],
               pch=21, col=col[wh], bg=bg[wh], ...)
        points(y.sub[!wh,1], y.sub[!wh,2],
               pch=21, col=col[!wh], bg=bg[!wh], ...)

        if(nrow(ynog) > 0) { # phenotype no genotype
            theinf <- ginf[rownames(ynog)]

            infcol <- rep("", length(nrow(ynog)))
            theinf[is.na(theinf)] <- 4
            for(i in seq(along=gnames))
                infcol[theinf==i] <- thecolors[i]
            obscol <- rep("black", length(infcol))

            col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
            bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

            points(ynog[,1], ynog[,2],
                   pch=21, col=col, bg=bg, ...)
        }
    }
    else { # pairs plot

        if(nrow(ynog) > 0) gnames <- c(gnames, "NA")

        # point colors
        obscol <- rep("", length(gobs.sub))
        gobs.sub[is.na(gobs.sub)] <- 4
        for(i in seq(along=gnames))
            obscol[gobs.sub==i] <- thecolors[i]
        infcol <- rep("", length(ginf.sub))
        ginf.sub[is.na(ginf.sub)] <- 4
        for(i in seq(along=gnames))
            infcol[ginf.sub==i] <- thecolors[i]

        if(nrow(ynog) > 0) { # phenotype no genotype
            theinf <- ginf[rownames(ynog)]

            infadd <- rep("", length(nrow(ynog)))
            theinf[is.na(theinf)] <- 4
            for(i in seq(along=gnames))
                infadd[theinf==i] <- thecolors[i]

            y.sub <- rbind(y.sub, ynog)
            obscol <- c(obscol, rep("black", nrow(ynog)))
            infcol <- c(infcol, infadd)
        }

        col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
        bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

        o <- order(as.numeric((col!=bg)))
        y.sub <- y.sub[o,,drop=FALSE]
        col <- col[o]
        bg <- bg[o]

        par(oma=c(0,0,1.5,0))
        pairs(y.sub, pch=21, col=col, bg=bg, ...)
        mtext(side=3, outer=TRUE, eqtlnam)
        wh <- (col == bg)
    }
}
