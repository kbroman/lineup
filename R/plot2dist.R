## plot2dist.R
## Karl W Broman

#' Plot two sets of inter-individual distances against one another
#'
#' Plot two sets of inter-individual distances against one another, colored by
#' self and non-self distances.
#'
#'
#' @param d1 Output of \code{\link{distee}}.
#' @param d2 Output of \code{\link{distee}}.
#' @param hirow Names of rows to highlight in green.
#' @param hicol Names of columns to highlight in orange.
#' @param xlab X-axis label (optional)
#' @param ylab Y-axis label (optional)
#' @param smoothScatter If TRUE, plot non-self distances with
#' \code{\link[graphics]{smoothScatter}}; if FALSE, use
#' \code{\link[graphics]{plot}}.
#' @param colself Color to use for the self-self points.  If NULL, these aren't
#' plotted.
#' @param colnonself Color to use for the non-self points.  If NULL, these
#' aren't plotted.
#' @param colhirow Color to use for the \code{hirow} points.  If NULL, these
#' aren't plotted.
#' @param colhicol Color to use for the \code{hicol} points.  If NULL, these
#' aren't plotted.
#' @param \dots Passed to \code{\link[graphics]{plot}} and
#' \code{\link[graphics]{points}}.
#' @return None.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{pulldiag}}, \code{\link{distee}},
#' \code{\link{summary.lineupdist}}
#' @keywords graphics
#' @examples
#'
#' \dontrun{
#' # simulate MVN, 100 individuals, 40 measurements (of which 20 are just noise)
#' V <- matrix(0.3, ncol=20, nrow=20) + diag(rep(0.5, 20))
#' D <- chol(V)
#' z <- matrix(rnorm(20*100), ncol=20) %*% D
#'
#' # create two data matrices as z + noise
#' x <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))
#' y <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))
#'
#' # permute some rows
#' x[51:53,] <- x[c(52,53,51),]
#' y[41:42,] <- y[42:41,]
#'
#' # add column and row names
#' dimnames(x) <- dimnames(y) <- list(paste("ind", 1:100, sep=""),
#'                                    paste("gene", 1:40, sep=""))
#'
#' # calculate correlations between cols of x and cols of y
#' thecor <- corbetw2mat(x, y)
#'
#' # subset x and y, taking only columns with corr > 0.75
#' xs <- x[,thecor > 0.8]
#' ys <- y[,thecor > 0.8]
#'
#' # calculate distance (using "RMS difference" as a measure)
#' d1 <- distee(xs, ys, d.method="rmsd", labels=c("x","y"))
#'
#' # calculate distance (using "correlation" as a measure...really similarity)
#' d2 <- distee(xs, ys, d.method="cor", labels=c("x", "y"))
#'
#' # pull out the smallest 8 self-self correlations
#' sort(pulldiag(d2))[1:8]
#'
#' # summary of results
#' summary(d1)
#' summary(d2)
#'
#' # order to put matches together
#' summary(d2, reorder="alignmatches")
#'
#' # plot histograms of RMS distances
#' plot(d1)
#'
#' # plot histograms of correlations
#' plot(d2)
#'
#' # plot distances against one another
#' plot2dist(d1, d2)
#' }
#'
#' @importFrom graphics plot points
#' @importFrom grDevices colorRampPalette
#' @export
plot2dist <-
    function(d1, d2, hirow, hicol, xlab, ylab, smoothScatter=FALSE,
             colself="black", colnonself="gray", colhirow="green", colhicol="orange", ...)
{
    xmis <- ymix <- FALSE
    if(missing(xlab)) {
        xmis <- TRUE
        meth <- attr(d1, "d.method")
        if(is.null(meth)) xlab <- "d1"
        else if(meth=="cor") xlab <- "Correlation"
        else if(meth=="rmsd") xlab <- "RMS difference"
        else if(meth=="prop.mismatch") xlab <- "Proportion mismatches"
        else xlab <- "d1"
    }
    if(missing(ylab)) {
        ymis <- TRUE
        meth <- attr(d2, "d.method")
        if(is.null(meth)) ylab <- "d2"
        else if(meth=="cor") ylab <- "Correlation"
        else if(meth=="rmsd") ylab <- "RMS difference"
        else if(meth=="prop.mismatch") ylab <- "Proportion mismatches"
        else ylab <- "d2"
    }
    if(xlab == ylab && (xmis || ymis)) {
        xlab <- "d1"
        ylab <- "d2"
    }

    if(any(dim(d1) != dim(d2)) || any(rownames(d1) != rownames(d2)) ||
       any(colnames(d1) != colnames(d2))) {
        # pull out just the common rows and columns
        rn1 <- rownames(d1)
        rn2 <- rownames(d2)
        cn1 <- colnames(d1)
        cn2 <- colnames(d2)

        # line up rows
        m1 <- match(rn1, rn2)
        m2 <- match(rn2, rn1)
        if(any(is.na(m1))) {
            d1 <- d1[!is.na(m1),,drop=FALSE]
            rn1 <- rn1[!is.na(m1)]
        }
        if(any(is.na(m2))) {
            d2 <- d2[!is.na(m2),,drop=FALSE]
            rn2 <- rn2[!is.na(m1)]
        }
        d2 <- d2[match(rn1, rn2),,drop=FALSE]

        # line up columns
        m1 <- match(cn1, cn2)
        m2 <- match(cn2, cn1)
        if(any(is.na(m1))) {
            d1 <- d1[,!is.na(m1),drop=FALSE]
            cn1 <- cn1[!is.na(m1)]
        }
        if(any(is.na(m2))) {
            d2 <- d2[,!is.na(m2),drop=FALSE]
            cn2 <- cn2[!is.na(m1)]
        }
        d2 <- d2[,match(cn1, cn2),drop=FALSE]
    }

    rn <- rownames(d1)
    cn <- colnames(d1)

    m <- match(rn, cn)
    self <- matrix(ncol=2, nrow=sum(!is.na(m)))
    wh <- which(!is.na(m))
    m <- m[!is.na(m)]
    xl <- range(d1, na.rm=TRUE)
    yl <- range(d2, na.rm=TRUE)
    for(i in seq(along=wh)) {
        self[i,] <- c(d1[wh[i],m[i]], d2[wh[i],m[i]])
        d1[wh[i],m[i]] <- d2[wh[i],m[i]] <- NA
    }
    if(!missing(hirow)) {
        hirowd1 <- d1[hirow,]
        hirowd2 <- d2[hirow,]
        d1[hirow,] <- d2[hirow,] <- NA
    }
    if(!missing(hicol)) {
        hicold1 <- d1[,hicol]
        hicold2 <- d2[,hicol]
        d1[,hicol] <- d2[,hicol] <- NA
    }
    if(is.null(colnonself))
        plot(0,0,type="n", xlab=xlab, ylab=ylab, xlim=xl, ylim=yl, ...)
    else {
        if(smoothScatter)
            graphics::smoothScatter(d1, d2, xlab=xlab, ylab=ylab, xlim=xl, ylim=yl,
                                    colramp=colorRampPalette(c("white","blue")))
        else {
            plot(unclass(d1), unclass(d2), xlab=xlab, ylab=ylab, xlim=xl, ylim=yl, col=colnonself, ...)
        }
    }
    if(!missing(hirow) && !is.null(colhirow)) points(hirowd1, hirowd2, col=colhirow, ...)
    if(!missing(hicol) && !is.null(colhicol)) points(hicold1, hicold2, col=colhicol, ...)
    if(!is.null(colself)) points(self, col=colself, pch=16, ...)
}
