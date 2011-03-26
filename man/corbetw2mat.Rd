\name{corbetw2mat}
\alias{corbetw2mat}

\title{Calculate correlations between columns of two matrices}

\description{
  For matrices x and y, calculate the correlation between columns of x
  and columns of y.
}

\usage{
corbetw2mat(x, y, paired=TRUE, scaled=FALSE, bestpairs=FALSE, corthresh=0.9)
}

\arguments{
  \item{x}{A numeric matrix.}
  \item{y}{A numeric matrix with the same number of rows as \code{x}.}
  \item{paired}{If TRUE, calculate correlation between each column of \code{x}
  and the corresponding column of \code{y}.  \code{x} and \code{y} must
  have the same number of columns.}
  \item{scaled}{If TRUE, assume that the columns of \code{x} and
    \code{y} each have mean 0 and SD 1.}
  \item{bestpairs}{If \code{paired=TRUE}, this is ignored.  Otherwise,
  if \code{bestpairs=TRUE}, we find all pairs of columns (with one in
  \code{x} and one in \code{y}) that have correlation \eqn{\ge}
  \code{corthresh}.  If \code{bestpairs=FALSE}, we consider each column
  of \code{x} and find the column of \code{y} with the highest
  correlation.}
  \item{corthresh}{Threshold on correlations if \code{paired=FALSE} and
  \code{bestpairs=TRUE}.}
}

\details{
  Missing values (\code{NA}) are ignored, and we calculate the
  correlation using all complete pairs, as in \code{\link[stats]{cor}}
  with \code{use="pairwise.complete.obs"}.

  If \code{scaled=TRUE}, we assume that the columns of \code{x} and
  \code{y} have already been standardized, to have mean 0 and SD 1.
  Otherwise, we use \code{\link{fscale}} to standardize the
  columns and then calculate correlations as
  \code{mean(x[,i]*y[,j])/(nrow(x)-1)}.
}

\value{
  If \code{paired=TRUE}, the return value is a vector of correlations, between
  columns of \code{x} and the corresponding column of \code{y}.
  \code{x} and \code{y} must have the same number of columns.

  If \code{paired=FALSE} and \code{bestpairs=FALSE}, we return a matrix
  of size \code{ncol(x)} by \code{2}, with the \eqn{i}th row being the maximum 
  correlation between a column \eqn{i} of \code{x} and a column of
  \code{y}, and then the \code{y}-column index with that correlation.
  (In case of ties, we give the first one.)

  If \code{paired=FALSE} and \code{bestpairs=TRUE}, we return a matrix
  with three columns, containing all pairs of columns (with one in
  \code{x} and one in \code{y}) with correlation \eqn{\ge} \code{corthresh}.
  Each row corresponds to a column pair, and contains the correlation
  and then the \code{x}- and \code{y}-column indices.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

%\references{
%}

\examples{
# a variance matrix
V <- diag(rep(0.5, 5)) + 0.5
D <- chol(V)

# simulate two correlated matrices
x <- matrix(rnorm(100), ncol=5) %*% D 
y <- matrix(rnorm(100), ncol=5) %*% D + x

# create shuffled version of the second matrix
u <- sample(1:ncol(y))
z <- y[,u]

# correlations with paired columns
corbetw2mat(x, y)

# the same with y columns shuffled
corbetw2mat(x, z)

# for each column x, find column of y with max correlation
corbetw2mat(x, y, paired=FALSE)

# the same with y columns shuffled 
corbetw2mat(x, z, paired=FALSE)


# all pairs of columns with correlation >= 0.6
corbetw2mat(x, y, paired=FALSE, bestpairs=TRUE, corthresh=0.6)

# the same with y columns shuffled 
corbetw2mat(x, z, paired=FALSE, bestpairs=TRUE, corthresh=0.6)
}

\seealso{ \code{\link{distee}}, \code{\link{fscale}}, \code{\link{findCommonID}} }

\keyword{array}
\keyword{univar}
\keyword{multivariate}

