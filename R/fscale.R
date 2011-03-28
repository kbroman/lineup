######################################################################
#
# fscale.R
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
# Contains: fscale
#
######################################################################


# standardize the columns of a matrix so that they have mean 0 and SD 1
fscale <-
function(x)
{
  if(!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  y <- matrix(.C("R_fscale",
                 as.integer(n),
                 as.integer(p),
                 x=as.double(x),
                 PACKAGE="lineup",
                 NAOK=TRUE)$x, ncol=p)
  dimnames(y) <- dimnames(x)
  y
}

# end of fscale.R
