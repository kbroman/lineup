######################################################################
#
# subset.lineupdist.R
#
# copyright (c) 2011, Karl W Broman
# last modified Apr, 2011
# first written Apr, 2011
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
# Contains: subset.lineupdist, [.lineupdist
#
######################################################################

##############################
# subset.lineupdist
##############################

subset.lineupdist <-
function(x, rows, cols, ...)
{
  if(missing(cols) && !missing(rows) && length(rows)==prod(dim(x)))
    return(unclass(x)[rows])

  xnew <- unclass(x)[rows,cols]

  class(xnew) <- class(x)

  possible.attributes <- c("d.method", "labels", "compareWithin", "orig.selfd",
                           "badind", "obsg", "infg", "y", "denom", "linkwts", "genonames")
  for(i in possible.attributes[possible.attributes %in% names(attributes(x))])
    attr(xnew, i) <- attr(x, i)

  xnew
}
  
`[.lineupdist` <-
function(x, rows, cols)
subset(x, rows, cols)  

# end of subset.lineupdist.R
