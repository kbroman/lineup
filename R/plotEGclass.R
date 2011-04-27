######################################################################
#
# plotEGclass.R
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
# Contains: plotEGclass
#
######################################################################

######################################################################
# plot expression -> eQTL genotype classifier
######################################################################
plotEGclass <-
function(d, eqtl=1, outercol="inferred", innercol="observed", 
         thecolors=c("blue","green","red","orange"), ...)
{
  if(!("eg.lineupdist" %in% class(d)))
    stop("Input d must be as produced by disteg().")

  # inner and outer colors based on...
  choices <- c("observed", "inferred")
  if(!is.na(pm <- pmatch(outercol, choices)))
    outercol <- choices[pm]
  if(!is.na(pm <- pmatch(innercol, choices)))
    innercol <- choices[pm]

  ginf <- attr(d, "infg")
  gobs <- attr(d, "obsg")
  y <- attr(d, "y")
  gnames <- attr(d, "genonames")
  
  if(length(eqtl) > 1) {
    warning("eqtl should have length 1; using first element.")
    eqtl <- eqtl[1]
  }
  if(is.character(eqtl)) { # name of an eQTL
    eqtlnam <- eqtl
    eqtl <- match(eqtl, colnames(ginf))
    if(is.na(eqtl))
      stop("Can't find eQTL ", eqtl)
  }
  else {
    if(eqtl < 1 || eqtl > ncol(ginf))
      stop("eqtl must be between 1 and ", ncol(ginf))

    eqtlnam <- colnames(ginf)[eqtl]
  }

  # clean up eQTL name
  if(length(grep("@", eqtlnam)) == 1) {
    spl <- unlist(strsplit(eqtlnam, "@"))
    splnum <- unlist(strsplit(as.character(round(as.numeric(spl[2]),2)), "\\."))
    if(length(splnum)==1 || nchar(splnum[2]) == 0)
      splnum[2] <- "00"
    else if(nchar(splnum[2]) == 1)
      splnum[2] <- paste(splnum[2], "0", sep="")
    eqtlnam <- paste(spl[1], " @ ", splnum[1], ".", splnum[2], sep="")
  }

  ids <- findCommonID(rownames(gobs), rownames(ginf))

  ginf <- ginf[,eqtl]
  gobs <- gobs[,eqtl]
  y <- y[[eqtl]]
  hasy <- rowSums(is.na(y)) == 0

  gobs.sub <- gobs[ids$first]
  ginf.sub <- ginf[ids$second]
  y.sub <- y[ids$second,,drop=FALSE]

  ynog <- y[-ids$second,,drop=FALSE]
  ynog <- ynog[rowSums(is.na(ynog))==0,,drop=FALSE]

  if(ncol(y)==1) { # dot plot
    if(nrow(ynog) > 0) gnames <- c(gnames, "NA")

    xjit <- 0.2
    u <- runif(length(gobs.sub), -xjit, xjit)

    # point colors
    obscol <- rep("", length(gobs.sub))
    gobs.sub[is.na(gobs.sub)] <- 4
    for(i in seq(along=gnames))
      obscol[gobs.sub==i] <- thecolors[i]
    infcol <- rep("", length(ginf.sub))
    ginf.sub[is.na(ginf.sub)] <- 4
    for(i in seq(along=gnames))
      infcol[ginf.sub==i] <- thecolors[i]

    col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
    bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

    plot(gobs.sub+u, y.sub, xlab="Observed genotype", ylab=colnames(y),
         main=eqtlnam, xaxt="n", xlim=c(0.5, length(gnames)+0.5),
         pch=21, col=col, bg=bg, ylim=range(y, na.rm=TRUE), type="n", ...)
    axis(side=1, at=seq(along=gnames), gnames)

    # make sure the mismatches are on top
    wh <- !is.na(gobs.sub) & !is.na(ginf.sub) & gobs.sub==ginf.sub
    points((gobs.sub+u)[wh], y.sub[wh],
           pch=21, col=col[wh], bg=bg[wh], ...)
    points((gobs.sub+u)[!wh], y.sub[!wh],
           pch=21, col=col[!wh], bg=bg[!wh], ...)

    if(nrow(ynog) > 0) { # phenotype no genotype
      theinf <- ginf[rownames(ynog)]

      infcol <- rep("", length(nrow(ynog)))
      theinf[is.na(theinf)] <- 4
      for(i in seq(along=gnames))
        infcol[theinf==i] <- thecolors[i]
      obscol <- rep("black", length(infcol))

      col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
      bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

      points(runif(nrow(ynog), 4-xjit, 4+xjit), ynog,
             pch=21, col=col, bg=bg, ...)
    }

  }
  else if(ncol(y) == 2) { # scatter plot

    if(nrow(ynog) > 0) gnames <- c(gnames, "NA")

    # point colors
    obscol <- rep("", length(gobs.sub))
    gobs.sub[is.na(gobs.sub)] <- 4
    for(i in seq(along=gnames))
      obscol[gobs.sub==i] <- thecolors[i]
    infcol <- rep("", length(ginf.sub))
    ginf.sub[is.na(ginf.sub)] <- 4
    for(i in seq(along=gnames))
      infcol[ginf.sub==i] <- thecolors[i]

    col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
    bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

    plot(y.sub[,1], y.sub[,2], xlab=colnames(y)[1], ylab=colnames(y)[2],
         main=eqtlnam, pch=21, col=col, bg=bg, 
         xlim=range(y[,1], na.rm=TRUE), ylim=range(y[,2], na.rm=TRUE), type="n", ...)

    # make sure the mismatches are on top
    wh <- !is.na(gobs.sub) & !is.na(ginf.sub) & gobs.sub==ginf.sub
    points(y.sub[wh,1], y.sub[wh,2],
           pch=21, col=col[wh], bg=bg[wh], ...)
    points(y.sub[!wh,1], y.sub[!wh,2],
            pch=21, col=col[!wh], bg=bg[!wh], ...)

    if(nrow(ynog) > 0) { # phenotype no genotype
      theinf <- ginf[rownames(ynog)]

      infcol <- rep("", length(nrow(ynog)))
      theinf[is.na(theinf)] <- 4
      for(i in seq(along=gnames))
        infcol[theinf==i] <- thecolors[i]
      obscol <- rep("black", length(infcol))

      col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
      bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

      points(ynog[,1], ynog[,2],
             pch=21, col=col, bg=bg, ...)
    }
  }
  else { # pairs plot

    if(nrow(ynog) > 0) gnames <- c(gnames, "NA")

    # point colors
    obscol <- rep("", length(gobs.sub))
    gobs.sub[is.na(gobs.sub)] <- 4
    for(i in seq(along=gnames))
      obscol[gobs.sub==i] <- thecolors[i]
    infcol <- rep("", length(ginf.sub))
    ginf.sub[is.na(ginf.sub)] <- 4
    for(i in seq(along=gnames))
      infcol[ginf.sub==i] <- thecolors[i]

    if(nrow(ynog) > 0) { # phenotype no genotype
      theinf <- ginf[rownames(ynog)]

      infadd <- rep("", length(nrow(ynog)))
      theinf[is.na(theinf)] <- 4
      for(i in seq(along=gnames))
        infadd[theinf==i] <- thecolors[i]

      y.sub <- rbind(y.sub, ynog)
      obscol <- c(obscol, rep("black", nrow(ynog)))
      infcol <- c(infcol, infadd)
    }

    col <- switch(outercol, "observed"=obscol, "inferred"=infcol, rep(outercol, length(obscol)))
    bg <- switch(innercol, "observed"=obscol, "inferred"=infcol, rep(innercol, length(obscol)))

    o <- order(as.numeric((col!=bg)))
    y.sub <- y.sub[o,,drop=FALSE]
    col <- col[o]
    bg <- bg[o]

    par(oma=c(0,0,1.5,0))
    pairs(y.sub, pch=21, col=col, bg=bg, ...)
    mtext(side=3, outer=TRUE, eqtlnam)
    wh <- (col == bg)
  }
}

# end of plotEGclass.R
