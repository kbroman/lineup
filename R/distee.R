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
# e1, e2 = the two data sets (individuals x genes, or if transpose=TRUE,
#          genes x individuals)
#
# cor.threshold = threshold on between-gene correlation; use genes 
#                 with between-condition correlation above this to
#                 calculate the distance
#
# n.col = if cor.threshold not given, instead provide the number of
#         genes to use (specified by this argument).  This number will
#         be chosen, in decreasing order of correlation
#
# d.method = calculate inter-individual distance as RMS difference or
#            as correlation
#
# transpose = if TRUE, e1 and e2 are input as genes x individuals and
#             need to be transposed
#
# labels = labels to attach to the two data sets
#
# weightByCorr = if TRUE, in the calculate of inter-individual
#                distances, weight genes by their cross-tissue
#                correlation
#
# standardize = if TRUE, standardize the columns of the expression
#               matrices so they have mean 0 and SD 1
#
# verbose = if TRUE, give verbose output
# 
######################################################################

distee <-
function(e1, e2, cor.threshold, n.col, d.method=c("rmsd", "cor"),
         transpose=FALSE, labels=c("e1","e2"), weightByCorr=FALSE,
         standardize=FALSE, verbose=TRUE)
{
  if(!missing(e2) && missing(cor.threshold) && missing(n.col))
    stop("Give either cor.threshold or n.col")
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
  if(d.method=="rmsd") {
    d.func <-
      function(a,b,weights=NULL) {
        if(is.null(weights)) return(sqrt(mean((a-b)^2, na.rm=TRUE)))
        wh <- !is.na(a) & !is.na(b)
        a <- a[wh]
        b <- b[wh]
        weights <- weights[wh]
        sqrt(sum((a-b)^2*weights)/sum(weights))
      }
  }
  else {
    d.func <-
      function(a,b,weights=NULL) {
        if(is.null(weights)) return(cor(a,b, use="complete"))
        wh <- !is.na(a) & !is.na(b)
        cov.wt(cbind(a[wh],b[wh]), wt=weights[wh], cor=TRUE)$cor[1,2]
      }
  }
  if(missing(e2)) { # need to have pre-selected genes
    compareWithin <- TRUE
    cor.threshold <- -2
    if(transpose) e1 <- t(e1)
    e2 <- e1
  }
  else {
    compareWithin <- FALSE
    
    if(transpose) {
      e1 <- t(e1)
      e2 <- t(e2)
    }
  }

  if(compareWithin) weightByCorr <- FALSE 

  if(any(colnames(e1) != colnames(e2))) {
    cn1 <- colnames(e1)
    cn2 <- colnames(e2)
    m1 <- match(cn1, cn2)
    if(any(is.na(m1))) {
      if(verbose) cat("Omitting", sum(is.na(m1)), "genes from e1\n")
      e1 <- e1[,!is.na(m1),drop=FALSE]
    }
    if(any(is.na(m2))) {
      if(verbose) cat("Omitting", sum(is.na(m2)), "genes from e2\n")
      e2 <- e2[,!is.na(m2),drop=FALSE]
    }
    m1 <- match(cn1[!is.na(m1)], cn2[!is.na(m2)])
    e2 <- e2[,m1]
    if(ncol(e1) < 3)
      stop("Need at least 3 genes in common (and should have many more).\n")
  }
  if(!missing(n.col) && n.col > ncol(e1)) {
    n.col <- ncol(e1)
    if(verbose) cat("Changing n.col to", n.col, "(all genes in common)\n")
  }
    
  if(standardize) {
    e1 <- scale(e1)
    if(!compareWithin) e2 <- scale(e2)
  }

  o1 <- e1 
  if(compareWithin) e2 <- e1
  o2 <- e2

  # first match cols and rows
  m1 <- match(rownames(e1), rownames(e2))
  m2 <- match(rownames(e2), rownames(e1))
  if(any(is.na(m1))) e1 <- e1[!is.na(m1),]
  if(any(is.na(m2))) e2 <- e2[!is.na(m2),]

  m1 <- match(rownames(e1), rownames(e2))

  e2 <- e2[m1,]

  if(any(rownames(e1) != rownames(e2)))
    stop("Mismatch in individual names")

  if(verbose) {
    cat("No. ind'ls = ", nrow(e1), "\n")
    cat("No. genes  = ", ncol(e1), "\n")
  }

  if((!missing(cor.threshold) && cor.threshold <= -1) ||
     (!missing(n.col) && n.col>=ncol(e1))) {
    keepAll <- TRUE
    if(verbose) cat("Retained all transcripts\n")
  }
  else keepAll <- FALSE

  if(!keepAll || weightByCorr) {
    thecor <- corbetw2mat(e1, e2, paired=TRUE,
                          scaled=standardize)

    if(missing(cor.threshold)) 
      cor.threshold <- sort(thecor, decreasing=TRUE)[n.col]

    if(!keepAll) {
      o1 <- o1[,thecor >= cor.threshold]
      o2 <- o2[,thecor >= cor.threshold]
    }
    if(!keepAll && weightByCorr) 
      thecor <- thecor[thecor >= cor.threshold]
    
    if(verbose) cat("Retained", ncol(o1), "transcripts\n")
  }
  if(!weightByCorr) thecor <- NULL

  d <- matrix(nrow=nrow(o1), ncol=nrow(o2))
  dimnames(d) <- list(rownames(o1), rownames(o2))
  if(compareWithin) {
    for(i in 1:(nrow(o1)-1))
      for(j in (i+1):nrow(o2)) 
        d[i,j] <- d[j,i] <- d.func(o1[i,], o2[j,], thecor)
  }
  else {
    for(i in 1:nrow(o1))
      for(j in 1:nrow(o2)) 
        d[i,j] <- d.func(o1[i,], o2[j,], thecor)
  }

  class(d) <- c("ee.lineupdist", "lineupdist")
  attr(d, "d.method") <- d.method
  attr(d, "retained") <- colnames(o1)
  attr(d, "labels") <- labels
  attr(d, "compareWithin") <- compareWithin
  d
}

# end of distee.R
