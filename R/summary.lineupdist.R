## summary.lineupdist.R
## Karl W Broman

# summary.lineupdist
#
#' Summarize inter-individual distances
#'
#' Summarize the results of \code{\link{distee}} or \code{\link{disteg}}, with
#' inter-individual distances between two sets of gene expression data.
#'
#'
#' @param object Output of \code{\link{distee}} or \code{\link{disteg}}.
#' @param cutoff (Optional) Cutoff on correlation/distance, with rows in the
#' results only being kept if the best distance/correlation is above this
#' cutoff or the self-self result is not missing and is above this cutoff.
#' @param dropmatches If TRUE, omit rows for which an individual's best match
#' is itself.
#' @param reorder If \code{"bydistance"}, reorder rows by increasing distance
#' (or decreasing correlation) to the best match and then by decreasing
#' distance (or decreasing correlation) to self; if \code{"alignmatches"},
#' group related errors together; if \code{"no"}, leave as is.
#' @param \dots Passed to \code{\link[base]{print.data.frame}}.
#' @return A list with two components: the distances summarized by row and the
#' distances summarized by column.
#'
#' For each individual, we calculate the minimum distance to others,
#' next-smallest distance, the self-self distance, the mean and SD of the
#' distances to others, and finally indicate the individual (or individuals)
#' that is closest.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{pulldiag}}, \code{\link{omitdiag}},
#' \code{\link{distee}}, \code{\link{disteg}}, \code{\link{plot2dist}},
#' \code{\link{plot.lineupdist}}
#' @keywords utilities
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
#' # summary of results, putting apparent matches together
#' summary(d1)
#' summary(d2)
#'
#' # order by correlations
#' summary(d2, reorder="bydistance")
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
#' @export
summary.lineupdist <-
    function(object, cutoff, dropmatches=TRUE, reorder=c("alignmatches", "bydistance", "no"), ...)
{
    d.method <- attr(object, "d.method")
    if(is.null(d.method)) d.method <- "rmsd"

    # cor: replace with negatives so sorting biggest to smallest
    if(d.method=="cor") object <- -object

    reorder <- match.arg(reorder)

    # comparison within tissue?
    compareWithin <- attr(object, "compareWithin")
    if(is.null(compareWithin)) compareWithin <- FALSE

    byrow <- data.frame(mind=apply(object, 1, min, na.rm=TRUE),
                        nextd=apply(object, 1, function(a) sort(a)[2]),
                        selfd=rep(NA, nrow(object)),
                        mean=rowMeans(object, na.rm=TRUE),
                        sd=apply(object, 1, sd, na.rm=TRUE),
                        best=apply(object, 1, function(a, b)
                        paste(b[!is.na(a) & a==min(a, na.rm=TRUE)], collapse=":"), colnames(object)))
    m <- match(rownames(object), colnames(object))
    for(i in which(!is.na(m)))
        byrow$selfd[i] <- object[i,m[i]]

    if(compareWithin) bycol <- byrow[-(1:nrow(byrow)),,drop=FALSE]
    else {
        bycol <- data.frame(mind=apply(object, 2, min, na.rm=TRUE),
                            nextd=apply(object, 2, function(a) sort(a)[2]),
                            selfd=rep(NA, ncol(object)),
                            mean=colMeans(object, na.rm=TRUE),
                            sd=apply(object, 2, sd, na.rm=TRUE),
                            best=apply(object, 2, function(a, b)
                            paste(b[!is.na(a) & a==min(a, na.rm=TRUE)], collapse=":"), rownames(object)))
        m <- match(colnames(object), rownames(object))
        for(i in which(!is.na(m)))
            bycol$selfd[i] <- object[m[i],i]
    }

    if(!compareWithin && dropmatches)
        res <- list(byrow=byrow[is.na(byrow$selfd) | byrow$selfd >= byrow$nextd,,drop=FALSE],
                    bycol=bycol[is.na(bycol$selfd) | bycol$selfd >= bycol$nextd,,drop=FALSE])
    else res <- list(byrow=byrow, bycol=bycol)

    if(reorder!="no") {
        if(compareWithin)
            res <- lapply(res, function(a) a[order(a$mind),])
        else
            res <- lapply(res, function(a) a[order(a$mind, -a$selfd),])
    }
    if(reorder=="alignmatches" && !compareWithin) {
        for(i in seq(along=res)) {
            if(nrow(res[[i]])<3) next
            col1 <- rownames(res[[i]])
            col2 <- as.character(res[[i]][,ncol(res[[i]])])
            n <- length(col1)
            theorder <- 1:n
            for(j in 1:(n-1)) {
                if(any(col2[theorder[j]] == col1[theorder[-(1:j)]])) {
                    wh <- which(col2[theorder[j]] == col1[theorder[-(1:j)]])

                    theorder[-(1:j)] <- theorder[c(wh+j, ((j+1):length(theorder))[-wh])]
                }
            }
            res[[i]] <- res[[i]][theorder,,drop=FALSE]
        }
    }

    # cor: return correlations to original scale and change colnames
    if(d.method=="cor") {
        res <- lapply(res, function(a) {a[,1:4] <- -a[,1:4]; a})
        for(i in 1:2) colnames(res[[i]])[1:3] <- c("maxc","nextc","selfc")
    }

    if(!missing(cutoff)) {
        if(d.method=="cor")
            res <- lapply(res, function(a,cut) a[a[,1]>=cut | (!is.na(a[,3]) & a[,3]>=cut),,drop=FALSE],cutoff)
        else
            res <- lapply(res, function(a,cut) a[a[,1]<=cut | (!is.na(a[,3]) & a[,3]<=cut),,drop=FALSE],cutoff)
    }

    class(res) <- "summary.lineupdist"
    attr(res, "labels") <- attr(object, "labels")
    attr(res, "d.method") <- attr(object, "d.method")
    attr(res, "retained") <- attr(object, "retained")
    attr(res, "compareWithin") <- compareWithin
    res
}

#' @export
print.summary.lineupdist <-
    function(x, ...)
{
    labels <- attr(x, "labels")
    if(is.null(labels))
        labels <- c("row", "col")

    cat("By ", labels[1], ":\n", sep="")
    print.data.frame(x$byrow, ...)

    compareWithin <- attr(x, "compareWithin")
    if(is.null(compareWithin)) compareWithin <- FALSE
    if(!compareWithin) {
        cat("\n")
        cat("By ", labels[2], ":\n", sep="")
        print.data.frame(x$bycol, ...)
    }
}

#' @export
print.lineupdist <-
    function(x, ...)
{
    possible.attributes <- c("d.method", "labels", "compareWithin", "orig.selfd",
                             "badind", "obsg", "infg", "y", "denom", "linkwts", "genonames")
    for(i in possible.attributes[possible.attributes %in% names(attributes(x))])
        attr(x, i) <- NULL

    print(unclass(x))
}
