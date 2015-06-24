## combinedist.R
## Karl W Broman

# combinedist
#
#' Combine distance matrices into a single such
#'
#' Combine multiple distance matrices into a single distance matrix providing
#' an overall summary
#'
#' The row and column names of the input distance matrices define the
#' individual IDs.
#'
#' If the input distance matrices all have an attribute \code{"denom"} (for
#' denominator) and \code{method="mean"}, we use a weighted mean, weighted by
#' the denominators.  This could be used to calculate an overall proportion.
#'
#' @param \dots Set of distance matrices, as calculated by \code{\link{distee}}
#' or \code{\link{disteg}}.
#' @param method Indicates whether to summarize using the median or the mean.
#' @return A distance matrix, with class \code{"lineupdist"}.  The individual
#' IDs are in the row and column names.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{distee}}, \code{\link{disteg}},
#' \code{\link{summary.lineupdist}}
#' @keywords utilities
#' @examples
#'
#' # simulate MVN, 100 individuals, 40 measurements (of which 20 are just noise)
#' V <- matrix(0.3, ncol=20, nrow=20) + diag(rep(0.5, 20))
#' D <- chol(V)
#' z <- matrix(rnorm(20*100), ncol=20) %*% D
#'
#' # create three data matrices as z + noise
#' x <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))
#' y <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))
#' w <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))
#'
#' # permute some rows of x
#' x[51:53,] <- x[c(52,53,51),]
#'
#' # add column and row names
#' dimnames(x) <- dimnames(y) <- dimnames(w) <-
#'    list(paste("ind", 1:100, sep=""), paste("gene", 1:40, sep=""))
#'
#' # calculate correlations between cols of x and of the other two matrices
#' corxy <- corbetw2mat(x, y)
#' corxw <- corbetw2mat(x, w)
#'
#' # using columns with corr > 0.75,
#' # calculate distance (using "correlation" as a measure...really similarity)
#' dxy <- distee(x[,corxy>0.75], y[,corxy>0.75], d.method="cor", labels=c("x", "y"))
#' dxw <- distee(x[,corxw>0.75], w[,corxw>0.75], d.method="cor", labels=c("x", "w"))
#'
#' d <- combinedist(dxy, dxw)
#'
#' summary(d)
#'
#' @export
combinedist <-
    function(..., method=c("median", "mean"))
{
    v <- list(...)

    # input is already a list?
    if(length(v) == 1 && is.list(v[[1]])) v <- v[[1]]

    if(!all(sapply(v, function(a) "lineupdist" %in% class(a))))
        stop("Input distance matrices must each be of class \"lineupdist\".")

    if(length(unique(sapply(v, function(a) class(a)[1]))) > 1)
        stop("Need all of the distance matrices to be the same type.")

    rn <- unique(unlist(lapply(v, rownames)))
    cn <- unique(unlist(lapply(v, colnames)))

    method <- match.arg(method)

    # combine into one big matrix
    d <- array(dim=c(length(rn), length(cn), length(v)))
    dimnames(d) <- list(rn, cn, names(v))
    for(i in seq(along=v))
        d[rownames(v[[i]]),colnames(v[[i]]),i] <- v[[i]]

    if(method=="mean" && all(sapply(v, function(a) !is.null(attr(a, "denom"))))) {
        use.denom <- TRUE
        denom <- array(dim=c(length(rn), length(cn), length(v)))
        dimnames(denom) <- list(rn, cn, names(v))
        for(i in seq(along=v))
            denom[rownames(v[[i]]), colnames(v[[i]]), i] <- attr(v[[i]], "denom")
    }
    else use.denom <- FALSE

    # summarize
    if(method=="median")
        ds <- apply(d, 1:2, median, na.rm=TRUE)
    else if(use.denom) {
        denom.sum <- apply(denom, 1:2, sum, na.rm=TRUE)
        ds <- apply(d*denom, 1:2, sum, na.rm=TRUE)/denom.sum
        attr(ds, "denom") <- denom.sum
    }
    else
        ds <- apply(d, 1:2, mean, na.rm=TRUE)

    class(ds) <- class(v[[1]])

    possible.attributes <- c("d.method", "compareWithin")
    for(i in possible.attributes[possible.attributes %in% names(attributes(v[[1]]))])
        attr(ds, i) <- attr(v[[1]], i)

    ds
}
