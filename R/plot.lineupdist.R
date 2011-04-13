######################################################################
#
# plot.lineupdist.R
#
# copyright (c) 2011, Karl W Broman
# last modified Apr, 2011
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
# Contains: plot.lineupdist
#
######################################################################

plot.lineupdist <-
function(x, breaks, add.rug=TRUE, what=c("both", "ss", "sn"), ...)
{
  what <- match.arg(what)
  if(what=="ns") what <- "sn"

  di <- pulldiag(x)
  if(what=="both") ra <- range(x, na.rm=TRUE)
  else if(what=="ss") ra <- range(di, na.rm=TRUE)
  else ra <- range(omitdiag(x), na.rm=TRUE)

  if(diff(ra)==0) ra <- ra+c(-0.001, 0.001)

  if(missing(breaks)) breaks <- seq(ra[1], ra[2], len=sqrt(prod(dim(x))))
  if(length(breaks)==1) breaks <- seq(ra[1], ra[2], len=breaks)
  d.method <- switch(attr(x, "d.method"), "cor"="correlation", "rmsd"="RMS distance")
  main <- paste(c("Self-self", "Self-nonself"), switch(attr(x, "d.method"), "cor"="correlation", "rmsd"="distance"))

  if(what != "sn" && !all(is.na(di))) { # some self-self distances

    if(what == "both") {
      old.mfrow <- par("mfrow")
      old.las <- par("las")
      on.exit(par(mfrow=old.mfrow, las=old.las))
      par(mfrow=c(2,1), las=1)
    }

    hist(di, breaks=breaks, xlab=d.method, main=main[1])
    if(add.rug) rug(di)
  }
  
  if(what != "ss") {
    x <- omitdiag(x)
    hist(x, breaks=breaks, xlab=d.method, main=main[2])
    if(add.rug) rug(x)
  }

  invisible()
}

# end of plot.lineupdist.R
