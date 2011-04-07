######################################################################
#
# omitdiag.R
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
# Contains: omitdiag
#
######################################################################

##############################
# replace "diagonal" (the self-self cases) from a distance matrix with NAs
##############################

omitdiag <-
function(d)
{
  rn <- rownames(d)
  cn <- colnames(d)
  m <- match(rn, cn)
  wh <- which(!is.na(m))
  m <- m[!is.na(m)]
  for(i in seq(along=wh)) 
    d[wh[i],m[i]] <- NA

  d
}
  
# end of omitdiag.R
