######################################################################
#
# plot.lineupdist.R
#
# copyright (c) 2011, Karl W Broman
# last modified Apr, 2011
# first written Mar, 2011
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/lineup package
# Contains: plot.lineupdist
#
######################################################################



#' Plot summary of inter-individual distances
#' 
#' Plot histograms of self-self and self-nonself distances from a distance
#' matrix calculated by \code{\link{distee}} or \code{\link{disteg}}.
#' 
#' We call \code{\link{pulldiag}} and \code{\link{omitdiag}} to get the
#' self-self and self-nonself distances.
#' 
#' If all of the self-self distances are missing, we plot just the self-nonself
#' distances.
#' 
#' @param x Output of \code{\link{distee}} or \code{\link{disteg}}.
#' @param breaks Optional vector of breaks, passed to
#' \code{\link[graphics]{hist}}, though if it is length 1, we interpret it as
#' the number of breaks and ensure that both histograms use the same set of
#' breaks.
#' @param add.rug If true, also include \code{\link[graphics]{rug}} below
#' histograms.
#' @param what Indicates whether to plot both self-self and self-nonself
#' distances (or correlations) or just one or the other.  (\code{"ss"}
#' indicates self-self and \code{"sn"} indicates self-nonself.)
#' @param \dots Ignored at this point.
#' @return None.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{pulldiag}}, \code{\link{distee}},
#' \code{\link{plot2dist}}
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
#' @importFrom graphics plot par hist
#' @method plot lineupdist
#' @export plot lineupdist
plot.lineupdist <-
function(x, breaks, add.rug=TRUE, what=c("both", "ss", "sn"), ...)
{
  what <- match.arg(what)
  if(what=="ns") what <- "sn"

  di <- pulldiag(x)
  if(what=="both") ra <- range(x, na.rm=TRUE)
  else if(what=="ss") ra <- range(di, na.rm=TRUE)
  else ra <- range(omitdiag(x), na.rm=TRUE)

  if(diff(ra)==0) ra <- ra+c(-0.001, 0.001)

  if(missing(breaks)) breaks <- seq(ra[1], ra[2], len=sqrt(prod(dim(x))))
  if(length(breaks)==1) breaks <- seq(ra[1], ra[2], len=breaks)
  d.method <- switch(attr(x, "d.method"), "cor"="correlation", "rmsd"="RMS distance")
  main <- paste(c("Self-self", "Self-nonself"), switch(attr(x, "d.method"), "cor"="correlation", "rmsd"="distance"))

  if(what != "sn" && !all(is.na(di))) { # some self-self distances

    if(what == "both") {
      old.mfrow <- par("mfrow")
      old.las <- par("las")
      on.exit(par(mfrow=old.mfrow, las=old.las))
      par(mfrow=c(2,1), las=1)
    }

    hist(di, breaks=breaks, xlab=d.method, main=main[1])
    if(add.rug) rug(di)
  }
  
  if(what != "ss") {
    x <- omitdiag(x)
    hist(x, breaks=breaks, xlab=d.method, main=main[2])
    if(add.rug) rug(x)
  }

  invisible()
}

# end of plot.lineupdist.R
