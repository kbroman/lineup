######################################################################
#
# distee.R
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
# Contains: distee
#
######################################################################

######################################################################
# distee: for calculating distance between pairs of individuals
#         with microarray-like data in two conditions (for example,
#         two tissues)
#
# e1, e2 = the two data sets (individuals x genes)
#
# d.method = calculate inter-individual distance as RMS difference or
#            as correlation
#
# labels = labels to attach to the two data sets
#
# weights = Optional vector of weights
#
# scaled = if TRUE, the columns of e1 and e2 are assumed to have
#          mean 0 and SD 1
#
# verbose = if TRUE, give verbose output
# 
######################################################################

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

  d <- matrix(nrow=nrow(e1), ncol=nrow(e2))
  dimnames(d) <- list(rownames(e1), rownames(e2))
  if(compareWithin) {
    for(i in 1:(nrow(o1)-1))
      for(j in (i+1):nrow(o2)) 
        d[i,j] <- d[j,i] <- d.func(e1[i,], e2[j,], thecor)
  }
  else {
    for(i in 1:nrow(e1))
      for(j in 1:nrow(e2)) 
        d[i,j] <- d.func(e1[i,], e2[j,], thecor)
  }

  class(d) <- c("ee.lineupdist", "lineupdist")
  attr(d, "d.method") <- d.method
  attr(d, "labels") <- labels
  attr(d, "compareWithin") <- compareWithin
  d
}

# end of distee.R
