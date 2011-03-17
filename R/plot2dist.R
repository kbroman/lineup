######################################################################
#
# plot2dist.R
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
# Contains: plot2dist
#
######################################################################

######################################################################
# plot two distances
######################################################################
plot2dist <-
function(d1, d2, hirow, hicol, xlab="RMS difference", ylab="correlation",
         smoothscatter=TRUE, ...)
{
  if(any(dim(d1) != dim(d2)) || any(rownames(d1) != rownames(d2)) ||
     any(colnames(d1) != colnames(d2))) {
    # pull out just the common rows and columns
    rn1 <- rownames(d1)
    rn2 <- rownames(d2)
    cn1 <- colnames(d1)
    cn2 <- colnames(d2)

    # line up rows
    m1 <- match(rn1, rn2)
    m2 <- match(rn2, rn1)
    if(any(is.na(m1))) {
      d1 <- d1[!is.na(m1),,drop=FALSE]
      rn1 <- rn1[!is.na(m1)]
    }
    if(any(is.na(m2))) {
      d2 <- d2[!is.na(m2),,drop=FALSE]
      rn2 <- rn2[!is.na(m1)]
    }
    d2 <- d2[match(rn1, rn2),,drop=FALSE]
    
    # line up columns
    m1 <- match(cn1, cn2)
    m2 <- match(cn2, cn1)
    if(any(is.na(m1))) {
      d1 <- d1[,!is.na(m1),drop=FALSE]
      cn1 <- cn1[!is.na(m1)]
    }
    if(any(is.na(m2))) {
      d2 <- d2[,!is.na(m2),drop=FALSE]
      cn2 <- cn2[!is.na(m1)]
    }
    d2 <- d2[,match(cn1, cn2),drop=FALSE]
  }

  rn <- rownames(d1)
  cn <- colnames(d1)

  m <- match(rn, cn)
  self <- matrix(ncol=2, nrow=sum(!is.na(m)))
  wh <- which(!is.na(m))
  m <- m[!is.na(m)]
  xl <- range(d1)
  yl <- range(d2)
  for(i in seq(along=wh)) {
    self[i,] <- c(d1[wh[i],m[i]], d2[wh[i],m[i]])
    d1[wh[i],m[i]] <- d2[wh[i],m[i]] <- NA
  }
  if(!missing(hirow)) {
    hirowd1 <- d1[hirow,]
    hirowd2 <- d2[hirow,]
    d1[hirow,] <- d2[hirow,] <- NA
  }
  if(!missing(hicol)) {
    hicold1 <- d1[,hicol]
    hicold2 <- d2[,hicol]
    d1[,hicol] <- d2[,hicol] <- NA
  }
  if(smoothscatter)
    smoothScatter(d1, d2, xlab=xlab, ylab=ylab, xlim=xl, ylim=yl,
                  colramp=colorRampPalette(c("white","blue")))
  else
    plot(d1, d2, xlab=xlab, ylab=ylab, xlim=xl, ylim=yl, col="gray", ...)
  if(!missing(hirow)) points(hirowd1, hirowd2, col="green", ...)
  if(!missing(hicol)) points(hicold1, hicold2, col="orange", ...)
  points(self, col="black", pch=16, ...)
}

# end of plot2dist.R
