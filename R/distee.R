## distee.R
## Karl W Broman

# distee
#
#' Calculate distance between two gene expression data sets
#'
#' Calculate a distance between all pairs of individuals for two gene
#' expression data sets
#'
#' We calculate the pairwise distance between all individuals (rows) in
#' \code{e1} and all individuals in \code{e2}.  This distance is either the RMS
#' difference (\code{d.method="rmsd"}) or the correlation
#' (\code{d.method="cor"}).
#'
#' @param e1 Numeric matrix of gene expression data, as individuals x genes.
#' The row and column names must contain individual and gene identifiers.
#' @param e2 (Optional) Like \code{e1}.  An appreciable number of individuals
#' and genes must be in common.
#' @param d.method Calculate inter-individual distance as RMS difference or as
#' correlation.
#' @param labels Two character strings, to use as labels for the two data
#' matrices in subsequent output.
#' @param verbose if TRUE, give verbose output.
#' @return A matrix with \code{nrow(e1)} rows and \code{nrow(e2)} columns,
#' containing the distances.  The individual IDs are in the row and column
#' names.  The matrix is assigned class \code{"lineupdist"}.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{pulldiag}}, \code{\link{omitdiag}},
#' \code{\link{summary.lineupdist}}, \code{\link{plot2dist}},
#' \code{\link{disteg}}, \code{\link{corbetw2mat}}
#' @keywords utilities
#' @examples
#'
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
#'
#' @useDynLib lineup
#' @export
distee <-
    function(e1, e2, d.method=c("rmsd", "cor"), labels=c("e1","e2"),
             verbose=TRUE)
{
    if(length(labels) != 2) {
        warning("labels should have length two; input ignored.")
        labels <- c("e1","e2")
    }
    if(is.null(colnames(e1)))
        stop("e1 is missing column names")
    if(is.null(rownames(e1)))
        stop("e1 is missing row names")
    if(!missing(e2)) {
        if(is.null(colnames(e2)))
            stop("e2 is missing column names")
        if(is.null(rownames(e2)))
            stop("e2 is missing row names")
    }

    d.method <- match.arg(d.method)

    if(missing(e2)) {
        e2 <- e1
        compareWithin <- TRUE
    }
    else {
        compareWithin <- FALSE

        # line up columns
        if(!compareWithin && ((ncol(e1) != ncol(e2)) ||
                              (colnames(e1) != colnames(e2)))) {
            cnmatch <- findCommonID(colnames(e1), colnames(e2))

            if(ncol(e1) != length(cnmatch$first)) {
                if(verbose) cat("Omitting", ncol(e1) - length(cnmatch$first), "genes from e1\n")
                e1 <- e1[,cnmatch$first,drop=FALSE]
            }
            if(ncol(e2) != length(cnmatch$second)) {
                if(verbose) cat("Omitting", ncol(e2) - length(cnmatch$second), "genes from e2\n")
                e2 <- e2[,cnmatch$second,drop=FALSE]
            }
        }
    }

    if(compareWithin) {
        if(d.method=="cor") {
            d <- cor(t(e1), use="pairwise.complete.obs")
            diag(d) <- NA
        }
        else
            d <- matrix(.C("R_rmsd",
                           as.integer(ncol(e1)),
                           as.integer(nrow(e1)),
                           as.double(t(e1)),
                           as.integer(nrow(e2)),
                           as.double(t(e2)),
                           d=as.double(rep(NA, nrow(e1)*nrow(e2))),
                           as.integer(1), # symmetric (e1==e2)
                           PACKAGE="lineup",
                           NAOK=TRUE)$d, ncol=nrow(e2))
    }
    else {
        if(d.method=="cor")
            d <- corbetw2mat(t(e1), t(e2), what="all")
        else
            d <- matrix(.C("R_rmsd",
                           as.integer(ncol(e1)),
                           as.integer(nrow(e1)),
                           as.double(t(e1)),
                           as.integer(nrow(e2)),
                           as.double(t(e2)),
                           d=as.double(rep(NA, nrow(e1)*nrow(e2))),
                           as.integer(0), # not symmetric (e1 != e2)
                           PACKAGE="lineup",
                           NAOK=TRUE)$d, ncol=nrow(e2))
    }

    dimnames(d) <- list(rownames(e1), rownames(e2))
    class(d) <- c("ee.lineupdist", "lineupdist")
    attr(d, "d.method") <- d.method
    attr(d, "labels") <- labels
    attr(d, "compareWithin") <- compareWithin
    d
}
