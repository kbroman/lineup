######################################################################
#
# disteg.R
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
# Contains: disteg
#
######################################################################

######################################################################
# disteg
#
# cross: A cross object 
#
# pheno: A matrix of expression phenotypes
#
# pmark: Pseudomarkers closest to each gene
#
# min.genoprob: Threshold on genotype probabilities; if maximum
#               probability is less than this, observed genotype taken
#               as NA.
# 
# k             Number of nearest neighbors to consider in forming a
#               k-nearest neighbor classifier
#
# min.classprob Minimum proportion of neighbors with a common class
#               to make a class prediction
#
# repeatKNN     If true, repeat k-nearest neighbor a second time,
#               after omitting individuals who seem to not be
#               self-self matches
#      
# min.selfd     Min distance from self (as proportion of mismatches
#               between observed and predicted eQTL genotypes) to be
#               excluded from the second round of k-nearest neighbor
#
# phenolabel    Label for expression phenotypes to place in the output
#               distance matrix
#
# weights       If not missing, use these as weights on the genes in
#               the calculation of inter-individual distances
#               (actually, since we're working with eQTL rather than
#               transcripts, we take the average weight for the
#               multiple transcripts that correspond to a give eQTL
#
# weightByLinkage If TRUE, weight the eQTL to account for their
#                 relative positions (for example, two tightly linked
#                 eQTL would each count about 1/2 of an isolated eQTL)
#
# map.function    Used if weightByLinkage is TRUE
#
# standardize = if TRUE, standardize the columns of the expression
#               matrices so they have mean 0 and SD 1
#
# verbose       If TRUE, print tracing information.
#
######################################################################

disteg <-
function(cross, pheno, pmark, min.genoprob=0.99,
         k=20, min.classprob=0.8, repeatKNN=TRUE,
         max.selfd=0.3, phenolabel="phenotype", weights,
         weightByLinkage=FALSE,
         map.function=c("haldane", "kosambi", "c-f", "morgan"),
         standardize=FALSE, verbose=TRUE)
{
  require(class)
  
  # individuals in common between two data sets
  theids <- findCommonID(cross, pheno)
  if(sum(theids$inBoth) < k)
    stop("You need at least ", k, " individuals in common between the data sets.")

  if(!missing(weights)) {
    if(ncol(pheno) != length(weights)) {
      weights <- rep(1, ncol(pheno))
      warning("weights should have length ncol(pheno), [", ncol(pheno),
              "]; ignored.")
    }
  }
  else
    weights <- rep(1, ncol(pheno))

  # make sure pheno and pmark line up
  m <- match(colnames(pheno), rownames(pmark))
  if(any(is.na(m))) pheno <- pheno[,!is.na(m),drop=FALSE]
  m <- match(rownames(pmark), colnames(pheno))
  if(any(is.na(m))) pmark <- pmark[!is.na(m),,drop=FALSE]
  m <- match(rownames(pmark), colnames(pheno))
  weights <- weights[m]
  pheno <- pheno[,m,drop=FALSE]
  
  # unique eQTL locations
  cpmark <- paste(pmark$chr, pmark$pmark, sep=":")
  upmark <- unique(cpmark)
  uweights <- rep(0, length(upmark))

  m <- match(upmark, cpmark)
  thechr <- pmark$chr[m]
  thepos <- pmark$pos[m]

  # construct weights to account for linkage
  if(weightByLinkage) {
    uchr <- unique(thechr)
    linkwts <- rep(1, length(thechr))

    map.function <- match.arg(map.function)
    mf <- switch(map.function,
                 "haldane" = mf.h,
                 "kosambi" = mf.k,
                 "c-f" = mf.cf,
                 "morgan" = mf.m)

    for(i in uchr) {
      d <- thepos[thechr==i]
      D <- matrix(1-2*mf(abs(outer(d, d, "-"))), ncol=length(d))
      linkwts[uchr==i] <- (length(d) + 1 - colSums(D))/length(d)
    }
  }



  # to contain observed and inferred eQTL genotypes
  obsg <- matrix(ncol=length(upmark), nrow=nind(cross))
  infg <- matrix(ncol=length(upmark), nrow=nrow(pheno))
  colnames(obsg) <- colnames(infg) <-
    apply(pmark[match(upmark, cpmark),c(1,3)], 1, paste, collapse="@")
  rownames(obsg) <- getid(cross)
  rownames(infg) <- rownames(pheno)
  
  # loop over eQTL; use k-nearest neighbor to classify
  if(verbose) cat("First pass through knn\n")
  for(i in seq(along=upmark)) {
    wh <- which(cpmark == upmark[i])
    uweights[i] <- mean(weights[cpmark==upmark[i]])
    
    y <- pheno[,wh,drop=FALSE]
    gp <- cross$geno[[pmark$chr[wh[1]]]]$prob[,pmark$pmark[wh[1]],,drop=FALSE]
    gi <- apply(gp, 1, function(a) which(a==max(a, na.rm=TRUE)))
    gmx <- apply(gp, 1, max, na.rm=TRUE)
    gi[gmx < min.genoprob] <- NA
    obsg[,i] <- gi

    if(standardize) {
      y <- scale(y)
      y[is.na(y)] <- 0
    }
    else {
      for(i in 1:ncol(y)) 
        y[is.na(y[,i]),i] <- mean(y[,i], na.rm=TRUE)
    }
    ysub <- y[theids[theids[,3],2],,drop=FALSE]
    gisub <- gi[theids[theids[,3],1]]
    keep <- !is.na(gisub) & apply(ysub, 1, function(a) !any(is.na(a) ))
    keep2 <- apply(y, 1, function(a) !any(is.na(a)))

    infg[!keep2,i] <- NA
    infg[keep2,i] <- knn(ysub[keep,,drop=FALSE], y[keep2,,drop=FALSE], gisub[keep],
                         k=k, l=ceiling(k*min.classprob))
  }
  if(weightByLinkage) uweights <- uweights * linkwts
  

  if(repeatKNN) {
    if(verbose) cat("Calculate self-self distances\n")
    # calculate self-self distances
    pd <- rep(NA, nrow(theids))
    names(pd) <- rownames(theids)
    for(i in rownames(theids)[theids[,3]]) 
      pd[i] <- mean(obsg[i,] != infg[i,], na.rm=TRUE)

    # bad individuals
    bad <- names(pd)[!is.na(pd) & pd>= max.selfd]
    subids <- theids[is.na(match(rownames(theids), bad)),,drop=FALSE]

    # repeat the k-nearest neighbor classification without the bad individuals
    if(verbose) cat("Second pass through knn\n")
    for(i in seq(along=upmark)) {
      wh <- which(cpmark == upmark[i])
      
      y <- pheno[,wh,drop=FALSE]
      gi <- obsg[,i]
  
      ysub <- y[subids[subids[,3],2],,drop=FALSE]
      gisub <- gi[subids[subids[,3],1]]
      keep <- !is.na(gisub) & apply(ysub, 1, function(a) !any(is.na(a) ))
      keep2 <- apply(y, 1, function(a) !any(is.na(a)))

      infg[!keep2,i] <- NA
      infg[keep2,i] <- knn(ysub[keep,,drop=FALSE], y[keep2,,drop=FALSE], gisub[keep],
                           k=k, l=ceiling(k*min.classprob))
    }
  }
  
  # calculate final distance
  if(verbose) cat("Calculate distance matrix\n")
  d <- matrix(nrow=nrow(obsg), ncol=nrow(infg))
  dimnames(d) <- list(rownames(obsg), rownames(infg))
  denom <- d
  for(i in 1:nrow(obsg)) {
    for(j in 1:nrow(infg)) {
      denom[i,j] <- sum((!is.na(obsg[i,]) & !is.na(infg[j,]))*uweights)
      d[i,j] <- sum((obsg[i,] != infg[j,])*uweights, na.rm=TRUE)/denom[i,j]
    }
  }

  attr(d, "d.method") <- "prop.mismatch"
  attr(d, "labels") <- c("genotype", phenolabel)
  attr(d, "retained") <- colnames(pheno)
  if(repeatKNN) attr(d, "orig.selfd") <- pd
  attr(d, "obsg") <- obsg
  attr(d, "infg") <- infg
  attr(d, "denom") <- denom
  names(uweights) <- upmark
  attr(d, "weights") <- uweights
  class(d) <- c("eg.lineupdist", "lineupdist")

  d
}

# end of disteg.R
