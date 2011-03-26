\name{findCommonID}
\alias{findCommonID}

\title{Find individuals in common between a cross and a phenotype matrix}

\description{
  Identify which individuals are in common between a QTL mapping data
  set and a matrix of phenotypes, series of genes.
}

\usage{
findCommonID(id1, id2)
}

\arguments{
 \item{id1}{A character vector of individual IDs.  This can also be a QTL cross
   object (see \code{\link[qtl]{read.cross}}), in which case
   \code{\link[qtl]{getid}} is used to grab individual IDs, or a matrix
   or data frame, in which case the rownames are taken to be IDs.}
 \item{id2}{Like \code{id1}, can be a character vector, a cross or a
 matrix/data frame.}
}

%\details{
%}

\value{
  A list with three components:
  
First, a data frame with rows corresponding to all individuals (across
the two sets of individual IDs) and three
columns: \code{indexInFirst} and \code{indexInSecond} contain numeric
indices to the locations of the individuals within \code{cross} and
\code{pheno}, and \code{inBoth} is a logical vector to indicate which
individuals appear in both crosses.  The row names are the
individual identifiers. 

The second and third components are vectors of indices in \code{id1} and
\code{id2}, respectively, indicating the paired locations of the
individuals that are in common.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

\examples{
id1 <- sample(LETTERS[1:5])
id2 <- LETTERS[3:8]
findCommonID(id1, id2)

x <- matrix(0, nrow=length(id2), ncol=3)
rownames(x) <- id2
findCommonID(id1, x)
}

\seealso{ \code{\link{calc.locallod}}, \code{\link{corbetw2mat}} }

\keyword{utilities}
