######################################################################
#
# calc.locallod.R
#
# copyright (c) 2011-2012, Karl W Broman
# last modified Apr, 2012
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
# Contains: calc.locallod
#
######################################################################

######################################################################
# calc.locallod
#
# cross: A cross object
#
# pheno: gene expression phenotypes (genes are rows)
#
# pmark: matrix containing chr and nearest pseudomarker (as from
#        find.gene.pseudomarker)
#
#
######################################################################

calc.locallod <-
function(cross, pheno, pmark, addcovar=NULL, intcovar=NULL, verbose=TRUE)
{
  require(qtl)

  if(any(pmark$chr == "X")) {
    warning("Dropping X chr loci; we can only handle autosomes for now.")
    pmark <- pmark[pmark$chr != "X",]
  }

  if(nind(cross) != nrow(pheno))
    stop("cross and pheno have incompatible numbers of individuals.")
  m <- match(colnames(pheno), rownames(pmark))
  if(any(is.na(m))) pheno <- pheno[,!is.na(m),drop=FALSE]
  m <- match(rownames(pmark), colnames(pheno))
  if(any(is.na(m))) pmark <- pmark[!is.na(m),,drop=FALSE]
  pheno <- pheno[,match(rownames(pmark), colnames(pheno)),drop=FALSE]

  cpmark <- paste(pmark$chr, pmark$pmark, sep=":")
  upmark <- unique(cpmark)

  lod <- rep(NA, ncol(pheno))
  temp <- cross
  n.ind <- nind(cross)

  # loop over unique pseudomarkers
  for(i in seq(along=upmark)) {
    if(verbose && i==round(i,-2)) cat(i, "of", length(upmark), "\n")
    wh <- which(cpmark == upmark[i])

    y <- pheno[,wh,drop=FALSE]
    gp <- cross$geno[[pmark$chr[wh[1]]]]$prob[,pmark$pmark[wh[1]],,drop=FALSE]

    # create dummy cross
    temp$pheno <- cbind(y, cross$pheno)
    temp$geno <- list("1"=list(data=cbind("m1"=rep(1, n.ind)), map=c("m1"=1),prob=gp))
    attr(temp$geno[[1]]$prob, "map") <- c("m1"=1)
    class(temp$geno[[1]]) <- "A"
    lod[wh] <- unlist(scanone(temp, method="hk", addcovar=addcovar, intcovar=intcovar,
                              pheno.col=1:ncol(y))[-(1:2)])
  }

  names(lod) <- colnames(pheno)
  lod
}


# end of calc.locallod.R
