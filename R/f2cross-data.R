#' Example experimental cross data
#'
#' Simulated experimental cross data with some sample mix-ups. The
#' only phenotype is an individual ID. There are 100 individuals
#' genotyped at 1000 markers on 19 autosomes.
#'
#' @docType data
#'
#' @usage data(f2cross)
#'
#' @format An object of class `"cross"`. See
#' [qtl::read.cross()] in the R/qtl package for details.
#'
#' @keywords datasets
#'
#' @seealso [expr1()], [expr2()], [genepos()], [pmap()]
#'
#' @examples
#' library(qtl)
#' data(f2cross)
#' summary(f2cross)
"f2cross"
