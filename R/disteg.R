######################################################################
#
# disteg.R
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
# classprob2drop If an individual is inferred to have a genotype
#                mismatch with classprob > this value, treat as an
#                outlier and drop from the analysis and then repeat
#                the KNN construction without it.
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
# weightByLinkage If TRUE, weight the eQTL to account for their
#                 relative positions (for example, two tightly linked
#                 eQTL would each count about 1/2 of an isolated eQTL)
#
# map.function    Used if weightByLinkage is TRUE
#
# verbose       If TRUE, print tracing information.
#
######################################################################

disteg <-
function(cross, pheno, pmark, min.genoprob=0.99,
         k=20, min.classprob=0.8, classprob2drop=1, repeatKNN=TRUE,
         max.selfd=0.3, phenolabel="phenotype", 
         weightByLinkage=FALSE,
         map.function=c("haldane", "kosambi", "c-f", "morgan"),
         verbose=TRUE)
{
  require(class)
  
  # individuals in common between two data sets
  theids <- findCommonID(cross, pheno)
  if(length(theids$first) < k)
    stop("You need at least ", k, " individuals in common between the data sets.")

  # make sure pheno and pmark line up
  m <- findCommonID(colnames(pheno), rownames(pmark))
  pheno <- pheno[,m$first,drop=FALSE]
  pmark <- pmark[m$second,,drop=FALSE]
  
  # drop X chromosome from cross
  chrtype <- sapply(cross$geno, class)
  if(any(chrtype=="X")) {
    warning("Dropping X chromosome")
    cross <- subset(cross, names(cross$geno)[chrtype=="A"])
  }
  crosschr <- names(cross$geno)

  # find transcript chr in cross
  m <- match(pmark$chr, crosschr)
  if(any(is.na(m))) {
    warning("Dropping ", sum(is.na(m)), " transcripts with unknown chromosome assignment.")
    pheno <- pheno[,!is.na(m), drop=FALSE]
    pmark <- pmark[!is.na(m),, drop=FALSE]
  }
  if(ncol(pheno) < 1)
    stop("Need at least one expression phenotype.")

  # unique eQTL locations
  cpmark <- paste(pmark$chr, pmark$pmark, sep=":")
  upmark <- unique(cpmark)

  m <- match(upmark, cpmark)
  thechr <- pmark$chr[m]
  thepos <- pmark$pos[m]

  # construct weights to account for linkage
  linkwts <- rep(1, length(thechr))
  if(weightByLinkage) {
    uchr <- unique(thechr)

    map.function <- match.arg(map.function)
    mf <- switch(map.function,
                 "haldane" = mf.h,
                 "kosambi" = mf.k,
                 "c-f" = mf.cf,
                 "morgan" = mf.m)

    for(i in uchr) {
      if(sum(thechr==i)==1) next # just one eQTL on this chromosome

      d <- thepos[thechr==i]
      D <- matrix(1-2*mf(abs(outer(d, d, "-"))), ncol=length(d))
      linkwts[thechr==i] <- (length(d) + 1 - colSums(D))/length(d)
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
  ysave <- vector("list", length(upmark))
  names(ysave) <- colnames(obsg)
  for(i in seq(along=upmark)) {
    wh <- which(cpmark == upmark[i])
    pmarkchr <- pmark$chr[wh[1]]
    pmarkpmark <- pmark$pmark[wh[1]]
    
    ysave[[i]] <- y <- pheno[,wh,drop=FALSE]
    gp <- cross$geno[[pmarkchr]]$prob[,pmarkpmark,,drop=FALSE]
    gi <- apply(gp, 1, function(a) which(a==max(a, na.rm=TRUE)))
    gmx <- apply(gp, 1, max, na.rm=TRUE)
    gi[gmx < min.genoprob] <- NA
    obsg[,i] <- gi

    ysub <- y[theids$second,,drop=FALSE]
    gisub <- gi[theids$first]
    keep <- !is.na(gisub) & (rowSums(is.na(ysub)) == 0) # have genotype and all phenotypes
    keep2 <- rowSums(is.na(y)) == 0 # have all phenotypes

    knnout <- knn(ysub[keep,,drop=FALSE], ysub[keep,,drop=FALSE], gisub[keep],
                  k=k, l=ceiling(k*min.classprob), prob=TRUE)
    pr <- attr(knnout, "prob")
    okeep <- keep
    keep[gisub[keep] != knnout & pr >= classprob2drop] <- FALSE
    if(verbose && sum(okeep) > sum(keep))
      cat(" -- Classifier ", i, ": dropping ", sum(okeep) - sum(keep), " outliers.\n", sep="")

    infg[keep2,i] <- knn(ysub[keep,,drop=FALSE], y[keep2,,drop=FALSE], gisub[keep],
                         k=k, l=ceiling(k*min.classprob))

  }


  if(repeatKNN) {
    if(verbose) cat("Calculate self-self distances\n")
    # calculate self-self distances
    pd <- rep(NA, nrow(theids$mat))
    names(pd) <- rownames(theids$mat)
    for(i in rownames(theids$mat)[theids$mat$inBoth])
      pd[i] <- mean(obsg[i,] != infg[i,], na.rm=TRUE)

    # bad individuals
    bad <- names(pd)[!is.na(pd) & pd>= max.selfd]

    # repeat the k-nearest neighbor classification without the bad individuals
    if(verbose) cat("Second pass through knn\n")
    for(i in seq(along=upmark)) {
      wh <- which(cpmark == upmark[i])
      
      y <- pheno[,wh,drop=FALSE]
      gi <- obsg[,i]
  
      ysub <- y[theids$second,,drop=FALSE]
      gisub <- gi[theids$first]
      keep <- !is.na(gisub) & (rowSums(is.na(ysub)) == 0) & is.na(match(rownames(ysub), bad))  
      keep2 <- rowSums(is.na(y)) == 0

      infg[!keep2,i] <- NA

      knnout <- knn(ysub[keep,,drop=FALSE], ysub[keep,,drop=FALSE], gisub[keep],
                    k=k, l=ceiling(k*min.classprob), prob=TRUE)
      pr <- attr(knnout, "prob")
      okeep <- keep
      keep[gisub[keep] != knnout & pr >= classprob2drop] <- FALSE
      if(verbose && sum(okeep) > sum(keep))
        cat(" -- Classifier ", i, ": dropping ", sum(okeep) - sum(keep), " outliers.\n", sep="")

      infg[keep2,i] <- knn(ysub[keep,,drop=FALSE], y[keep2,,drop=FALSE], gisub[keep],
                           k=k, l=ceiling(k*min.classprob))
    }
  }
  
  # calculate final distance
  if(verbose) cat("Calculate distance matrix\n")

  z <- .C("R_propmismatch",
          as.integer(ncol(obsg)),
          as.integer(nrow(obsg)),
          as.integer(t(obsg)),
          as.integer(nrow(infg)),
          as.integer(t(infg)),
          as.double(linkwts),
          prop=as.double(rep(NA, nrow(obsg)*nrow(infg))),
          denom=as.double(rep(NA, nrow(obsg)*nrow(infg))),
          PACKAGE="lineup",
          NAOK=TRUE)

  d <- matrix(z$prop, ncol=nrow(infg))
  denom <- matrix(z$denom, ncol=nrow(infg))
  dimnames(denom) <- dimnames(d) <- list(rownames(obsg), rownames(infg))

  attr(d, "d.method") <- "prop.mismatch"
  attr(d, "labels") <- c("genotype", phenolabel)
  if(repeatKNN) {
    attr(d, "orig.selfd") <- pd
    attr(d, "badind") <- bad
  }
  attr(d, "obsg") <- obsg
  attr(d, "infg") <- infg
  attr(d, "y") <- ysave
  attr(d, "denom") <- denom
  attr(d, "linkwts") <- linkwts
  attr(d, "genonames") <- getgenonames(class(cross)[1], "A", "simple",
                                       getsex(cross), attributes(cross))
  class(d) <- c("eg.lineupdist", "lineupdist")

  d
}

# end of disteg.R
