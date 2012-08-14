######################################################################
#
# combinedist.R
#
# copyright (c) 2012, Karl W Broman
# last modified Aug, 2012
# first written Aug, 2012
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
# Contains: combinedist
#
######################################################################

######################################################################
# combinedist: combine distance matrices into a single distance matrix
#
# ...: a set of distance matrices
#
# method = indicates whether to use median or mean
#
# if method = "mean" and attributes include "denom", we use a weighted
# mean, weighted by the denominator (e.g., to give overall proportions)
######################################################################

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

# end of combinedist.R
