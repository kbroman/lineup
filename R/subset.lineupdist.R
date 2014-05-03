## subset.lineupdist.R
## Karl W Broman

#' Subsetting distance matrix
#' 
#' Pull out a specified set of rows and columns from a distance matrix
#' calculated by \code{\link{distee}} or \code{\link{disteg}}.
#' 
#' 
#' @aliases subset.lineupdist [.lineupdist
#' @param x A distance matrix object as obtained from \code{\link{distee}} or
#' \code{\link{disteg}}.
#' @param rows Optional vector of selected rows.
#' @param cols Optional vector of selected columns.
#' @param \dots Ignored at this point.
#' @return The input distance matrix object, but with only the specified subset
#' of the data.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{disteg}}, \code{\link{distee}}, \code{\link{pulldiag}}
#' @keywords manip
#' @export
#' @rdname subset.lineupdist
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
  
#' @name [.lineupdist
#' @rdname subset.lineupdist
#' @export
`[.lineupdist` <-
function(x, rows, cols)
subset(x, rows, cols)  
