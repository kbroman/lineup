## util.R
## Karl W Broman

#' Installed version of R/lineup
#'
#' Print the version number of the currently installed version of R/lineup.
#'
#'
#' @return A character string with the version number of the currently
#' installed version of R/lineup.
#' @author Karl W Broman, \email{broman@@wisc.edu}
#' @keywords print
#' @examples
#'
#'   lineupversion()
#'
#' @export
lineupversion <-
    function()
{
    version <- unlist(utils::packageVersion("lineup"))

    if(length(version) == 3) {
        # make it like #.#-#
        return( paste(c(version, ".", "-")[c(1,4,2,5,3)], collapse="") )
    }

    paste(version, collapse=".")
}
