######################################################################
#
# corbetw2mat.R
#
# copyright (c) 2011, Karl W Broman
# last modified May, 2011
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
# Contains: corbetw2mat
#
######################################################################

######################################################################
# corbetw2mat
#
# x, y = two matrices with the same number of rows
#
# what = 
#  "paired":    We must have ncol(x) = ncol(y), and we then calculate 
#               the correlation for each column of x with the 
#               corresponding column of y
#
#  "bestright": We consider each column of x in turn, and find the
#               most column in y with the highest correlation,
#               returning a matrix with the max corr and the column 
#               index in y
#
#  "bestpairs": We find all pairs of columns (one in x, one in y)
#               with correlation >= corthreshold
#             
#  "all":       Calculate correlation betwen all pairs of columns
#               (one in x, one in y)
#
# corthresh: threshold on correlations to return; used only if
#            what="bestpairs"
#               
######################################################################

corbetw2mat <-
function(x, y, what=c("paired", "bestright", "bestpairs", "all"),
         corthresh=0.9)
{
  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.matrix(y)) y <- as.matrix(y)

  n <- nrow(x)
  if(nrow(y) != n)
    stop("nrow(x)=", n, ", which is not equal to nrow(y)=", nrow(y))

  px <- ncol(x)
  py <- ncol(y)
  what <- match.arg(what)
  
  if(is.null(colnames(x))) colnames(x) <- paste("V", 1:ncol(x), sep="")
  if(is.null(colnames(y))) colnames(y) <- paste("V", 1:ncol(y), sep="")

  if(what=="paired" && py != px)
    stop("what=\"paired\", but ncol(x)=", px, ", which is not equal to ncol(y)=", py)

  if(what=="paired") {
    res <- .C("R_corbetw2mat_paired",
              as.integer(n),
              as.integer(px),
              as.double(x),
              as.double(y),
              cor=as.double(rep(NA, px)),
              PACKAGE="lineup",
              NAOK=TRUE)$cor
    names(res) <- colnames(x)
  }

  else if(what=="bestright") {
    res <- .C("R_corbetw2mat_unpaired_lr",
              as.integer(n),
              as.integer(px),
              as.double(x),
              as.integer(py),
              as.double(y),
              cor=as.double(rep(NA, px)),
              index=as.integer(rep(NA, px)),
              PACKAGE="lineup",
              NAOK=TRUE)
    res <- data.frame(cor=res$cor, yindex=res$index)
    rownames(res) <- colnames(x)
    res <- cbind(res, ycol=colnames(y)[res[,2]])
  }
  else if(what=="bestpairs") {
    res <- .C("R_corbetw2mat_unpaired_best",
              as.integer(n),
              as.integer(px),
              as.double(x),
              as.integer(py),
              as.double(y),
              cor=as.double(rep(NA, px*py)),
              xindex=as.integer(rep(NA, px*py)),
              yindex=as.integer(rep(NA, px*py)),
              numpairs=as.integer(0),
              as.double(corthresh),
              PACKAGE="lineup",
              NAOK=TRUE)
    res <- data.frame(cor=res$cor[1:res$numpairs],
                 xindex=res$xindex[1:res$numpairs],
                 yindex=res$yindex[1:res$numpairs])
    res <- cbind(res, xcol=colnames(x)[res[,2]], ycol=colnames(y)[res[,3]])

  }
  else {
    res <- .C("R_corbetw2mat_unpaired_all",
              as.integer(n),
              as.integer(px),
              as.double(x),
              as.integer(py),
              as.double(y),
              cor=as.double(rep(NA, px*py)),
              PACKAGE="lineup",
              NAOK=TRUE)$cor
    res <- matrix(res, nrow=px, ncol=py)
    dimnames(res) <- list(colnames(x), colnames(y))
  }

  res
}

# end of corbetw2mat.R
