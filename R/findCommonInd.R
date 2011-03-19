######################################################################
#
# findCommonInd.R
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
# Contains: findCommonInd
#
######################################################################

######################################################################
# findCommonInd
#
# cross: A cross object 
#
# pheno: A matrix of phenotypes (row names = IDs)
#
# crossID: could be missing (in which case we use getid(cross))
#          if length 1, taken to be a phenotype name in cross that contains IDs
#          if length nind(cross), taken to be the cross IDs
#
######################################################################

findCommonInd <-
function(cross, pheno)
{
  crossID <- getid(cross)
  pheID <- rownames(pheno)
  
  m1 <- match(pheID, crossID)
  m2 <- match(crossID, pheID)

  totind <- length(crossID) + sum(is.na(m1))
  allID <- data.frame(indexInCross=rep(0, totind),
                      indexInPheno=rep(0, totind),
                      inBoth=rep(FALSE,totind))
  rownames(allID) <- c(crossID, pheID[is.na(m1)])

  allID[,1] <- match(rownames(allID), crossID)
  allID[,2] <- match(rownames(allID), pheID)
  allID[,3] <- !is.na(allID[,1]) & !is.na(allID[,2])
  allID
}

# end of findCommonInd.R
