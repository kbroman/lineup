\name{findCommonInd}
\alias{findCommonInd}

\title{Find individuals in common between a cross and a phenotype matrix}

\description{
  Identify which individuals are in common between a QTL mapping data
  set and a matrix of phenotypes, series of genes.
}

\usage{
findCommonInd(cross, pheno)
}

\arguments{
 \item{cross}{An object of class \code{"cross"} containing data for a
 QTL experiment.  See the help file for \code{\link[qtl]{read.cross}} in
 the R/qtl package (\url{http://www.rqtl.org}).  There must be a
 phenotype named \code{"id"} or \code{"ID"} that contains the
 individual identifiers.}
 \item{pheno}{A data frame of phenotypes, stored as individuals x
 phenotypes.  The row names must contain individual identifiers.}
}

%\details{
%}

\value{
A data frame with rows correspond to all individuals and three
columns: \code{indexInCross} and \code{indexInPheno} contain numeric
indices to the locations of the individuals within \code{cross} and
\code{pheno}, and \code{inBoth} is a logical vector to indicate which
individuals appear in both crosses.  The row names are the
individual identifiers. 
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

%\examples{
%}

\seealso{ \code{\link{calc.locallod}} }

\keyword{utilities}
