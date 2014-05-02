######################################################################
#
# pulldiag.R
#
# copyright (c) 2011-2012, Karl W Broman
# last modified Aug, 2012
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
# Contains: pulldiag
#
######################################################################

##############################
# pull out the "diagonal" from a distance matrix
#     (the self-self cases)
##############################



#' Pull out the diagonal from a distance matrix
#' 
#' Pull out the diagonal from a distance matrix calculated by
#' \code{\link{distee}} (that is, self-self distances).
#' 
#' We use the row and column names to identify which entries are self-self.
#' 
#' @param d A distance matrix calculated by \code{\link{distee}}.
#' @return A vector with the self-self distances.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{omitdiag}}, \code{\link{distee}}, \code{\link{disteg}},
#' \code{\link{summary.lineupdist}}, \code{\link{plot2dist}},
#' \code{\link{plot.lineupdist}}
#' @keywords array
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
#' @export pulldiag
pulldiag <-
function(d)
{
  ind <- findCommonID(rownames(d), colnames(d))
  diag(unclass(d)[ind$first,ind$second])
}
  
# end of pulldiag.R
