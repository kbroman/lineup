## pulldiag.R
## Karl W Broman

#' Pull out the diagonal from a distance matrix
#'
#' Pull out the diagonal from a distance matrix calculated by
#' [distee()] (that is, self-self distances).
#'
#' We use the row and column names to identify which entries are self-self.
#'
#' @param d A distance matrix calculated by [distee()].
#' @return A vector with the self-self distances.
#' @author Karl W Broman, \email{broman@@wisc.edu}
#' @seealso [omitdiag()], [distee()], [disteg()],
#' [summary.lineupdist()], [plot2dist()],
#' [plot.lineupdist()]
#' @keywords array
#' @examples
#' data(expr1, expr2)
#'
#' \dontshow{expr1 <- expr1[,1:500]
#' expr2 <- expr2[,1:500]}
#'
#' # distance as RMS difference
#' d <- distee(expr1, expr2)
#'
#' # pull out the self-self distances
#' d_selfself <- pulldiag(d)
#'
#' # samples with smallest self-self correlation
#' sort(d_selfself)[1:10]
#'
#' @export
pulldiag <-
    function(d)
{
    ind <- findCommonID(rownames(d), colnames(d))
    diag(unclass(d)[ind$first,ind$second])
}
