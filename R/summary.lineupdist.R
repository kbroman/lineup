######################################################################
#
# summary.lineupdist.R
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
# Contains: summary.lineupdist, print.summary.lineupdist
#
######################################################################

summary.lineupdist <-
function(object, dropmatches=TRUE, reorder=TRUE, ...)
{
  d.method <- attr(object, "d.method")
  if(is.null(d.method)) d.method <- "rmsd"

  # cor: replace with negatives so sorting biggest to smallest
  if(d.method=="cor") object <- -object

  bycol <- data.frame(mind=apply(object, 2, min),
                      nextd=apply(object, 2, function(a) sort(a)[2]),
                      selfd=rep(NA, ncol(object)),
                      mean=apply(object, 2, mean, na.rm=TRUE),
                      sd=apply(object, 2, sd, na.rm=TRUE),
                      best=apply(object, 2, function(a, b) paste(b[a==min(a)], collapse=":"), rownames(object)))
  m <- match(colnames(object), rownames(object))
  for(i in which(!is.na(m))) 
    bycol$selfd[i] <- object[m[i],i]

  compareWithin <- attr(object, "compareWithin")
  if(is.null(compareWithin)) compareWithin <- FALSE
  if(compareWithin) byrow <- bycol[-(1:nrow(bycol)),,drop=FALSE]
  else {
    byrow <- data.frame(mind=apply(object, 1, min),
                        nextd=apply(object, 1, function(a) sort(a)[2]),
                        selfd=rep(NA, nrow(object)),
                        mean=apply(object, 1, mean, na.rm=TRUE),
                        sd=apply(object, 1, sd, na.rm=TRUE),
                        best=apply(object, 1, function(a, b) paste(b[a==min(a)], collapse=":"), colnames(object)))
    m <- match(rownames(object), colnames(object))
    for(i in which(!is.na(m))) 
      byrow$selfd[i] <- object[i,m[i]]
  }

  if(dropmatches) 
    res <- list(byrow=byrow[is.na(byrow$selfd) | byrow$selfd >= byrow$nextd,,drop=FALSE],
                bycol=bycol[is.na(bycol$selfd) | bycol$selfd >= bycol$nextd,,drop=FALSE])
  else res <- list(byrow=byrow, bycol=bycol)
  
  if(reorder)
    res <- lapply(res, function(a) a[order(a$mind, -a$selfd),])
  
  # cor: return correlations to original scale and change colnames
  if(d.method=="cor") {
    res <- lapply(res, function(a) {a[,1:4] <- -a[,1:4]; a})
    for(i in 1:2) colnames(res[[i]])[1:3] <- c("maxc","nextc","selfc")
  }

  class(res) <- "summary.lineupdist"
  attr(res, "labels") <- attr(object, "labels")
  attr(res, "d.method") <- attr(object, "d.method")
  attr(res, "retained") <- attr(object, "retained")
  attr(res, "compareWithin") <- compareWithin
  res
}

print.summary.lineupdist <-
function(x, ...)
{
  labels <- attr(x, "labels")
  if(is.null(labels))
    labels <- c("row", "col")

  cat("By ", labels[1], ":\n", sep="")
  print.data.frame(x$byrow, ...)
  
  compareWithin <- attr(x, "compareWithin")
  if(is.null(compareWithin)) compareWithin <- FALSE
  if(!compareWithin) {
    cat("\n")
    cat("By ", labels[2], ":\n", sep="")
    print.data.frame(x$bycol, ...)
  }
}

# end of summary.lineupdist.R
