######################################################################
#
# pulldiag.R
#
# copyright (c) 2011-2012, Karl W Broman
# last modified Aug, 2012
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
# Contains: pulldiag
#
######################################################################

##############################
# pull out the "diagonal" from a distance matrix
#     (the self-self cases)
##############################

pulldiag <-
function(d)
{
  ind <- findCommonID(rownames(d), colnames(d))
  diag(unclass(d)[ind$first,ind$second])
}
  
# end of pulldiag.R
