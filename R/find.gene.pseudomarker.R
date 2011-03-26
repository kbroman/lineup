######################################################################
#
# find.gene.pseudomarker.R
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
# Contains: find.gene.pseudomarker
#
######################################################################

######################################################################
# find.gene.pseudomarker:
#         Function for identifying the closest pseudomarker to each
#         of a set of genes
#
# cross: A cross object (containing results of calc.genoprob or sim.geno)
#
# pmap:  Map object with physical locations of markers (positions in Mbp)
#
# geneloc: data.frame with chr + physical location (in Mbp) of genes
#          (gene names as row names)
#
# where: look in genotype probabilities or imputated genotypes?
#
######################################################################

find.gene.pseudomarker <-
function(cross, pmap, geneloc, where=c("prob", "draws"))
{
  where <- match.arg(where)
  if(!(where %in% names(cross$geno[[1]]))) 
    stop("You first need to run ", ifelse(where=="prob", "calc.genoprob", "sim.geno"), ".")
  
  require(qtl)

  cross <- replacemap(cross, pmap)
  res <- data.frame(chr=geneloc$chr,
                    pmark=find.pseudomarker(cross, geneloc$chr, geneloc$pos, where, addchr=FALSE),
                    stringsAsFactors=FALSE)

  rownames(res) <- rownames(geneloc)

  pmark <- res$pmark
  gr <- grep("^loc[0-9]+\\.*[0-9]*(\\.[0-9]+)*$", pmark)
  if(length(gr)>0) 
    pmark[gr] <- paste("c", res$chr[gr], ".", pmark[gr], sep="")
  upmark <- unique(pmark)
  thepos <- find.pseudomarkerpos(cross, upmark, where)
  res$pos <- thepos[match(pmark, rownames(thepos)),2]

  res <- cbind(res, dist.from.gene=(d <- geneloc$pos - res$pos))
  d <- abs(d)
  if(any(d > 2)) {
    ngap <- sum(d>2)
    maxd <- max(d)
    warning(ngap, " genes differ from pseudomarker pos by > 2 Mbp, with gaps as big as ", round(maxd, 1), " Mbp")
  }

  res
}

# end of find.gene.pseudomarker.R
