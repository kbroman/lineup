######################################################################
#
# corbetw2mat.R
#
# copyright (c) 2011, Karl W Broman
# last modified Mar, 2011
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
# paired: if TRUE, we must have ncol(x) = ncol(y), and we then 
#                  calculate the cor for each column of x with the 
#                  corresponding column of y
#
#         if FALSE...
#            if bestpairs=FALSE:
#                   we consider each column of x in turn, and find
#                   the most column in y with the highest correlation,
#                   returning a matrix with the max corr and the
#                   column index in y
#            if bestpairs=TRUE:
#                   we find all pairs of columns (one in x, one in y)
#                   with correlation >= corthreshold
#             
#
# scaled: if TRUE, we assume that the columns of each matrix has mean
#                  0 and SD 1, so we need only calculate
#                  sum(x[,i]*y[,i])/(n-1)
#
# bestpairs: only relevant if pairs=FALSE; see explanation above with
#            paired=FALSE
#            
# corthresh: threshold on correlations to return; used only if
#            pairs=FALSE and bestpairs=TRUE
#               
######################################################################

corbetw2mat <-
function(x, y, paired=TRUE, scaled=FALSE,
         bestpairs=FALSE, corthresh=0.9)
{
  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.matrix(y)) y <- as.matrix(y)

  n <- nrow(x)
  if(nrow(y) != n)
    stop("nrow(x)=", n, ", which is not equal to nrow(y)=", nrow(y))

  px <- ncol(x)
  py <- ncol(y)
  if(paired && py != px)
    stop("paired=TRUE, but ncol(x)=", px, ", which is not equal to ncol(y)=", py)

  if(paired & bestpairs)
    warning("bestpairs ignored when paired=TRUE")

  if(paired) 
    return(.C("R_corbetw2mat_paired",
              as.integer(n),
              as.integer(px),
              as.double(x),
              as.double(y),
              as.integer(scaled),
              cor=as.double(rep(NA, px)),
              PACKAGE="lineup",
              NAOK=TRUE)$cor)

  else {
    if(!bestpairs) {
      z <- .C("R_corbetw2mat_unpaired_lr",
              as.integer(n),
              as.integer(px),
              as.double(x),
              as.integer(py),
              as.double(y),
              as.integer(scaled),
              cor=as.double(rep(NA, px)),
              index=as.integer(rep(NA, px)),
              PACKAGE="lineup",
              NAOK=TRUE)
      return(cbind(cor=z$cor, yindex=z$index))
    }
    else {
      z <- .C("R_corbetw2mat_unpaired_best",
              as.integer(n),
              as.integer(px),
              as.double(x),
              as.integer(py),
              as.double(y),
              as.integer(scaled),
              cor=as.double(rep(NA, px*py)),
              xindex=as.integer(rep(NA, px*py)),
              yindex=as.integer(rep(NA, px*py)),
              numpairs=as.integer(0),
              as.double(corthresh),
              PACKAGE="lineup",
              NAOK=TRUE)
      return(cbind(cor=z$cor[1:z$numpairs],
                   xindex=z$xindex[1:z$numpairs],
                   yindex=z$yindex[1:z$numpairs]))
    }
  }
}

# end of corbetw2mat.R
