\name{corbetw2mat}
\alias{corbetw2mat}

\title{Calculate correlations between columns of two matrices}

\description{
  For matrices x and y, calculate the correlation between columns of x
  and columns of y.
}

\usage{
corbetw2mat(x, y, what=c("paired", "bestright", "bestpairs", "all"),
            scaled=FALSE, corthresh=0.9)
}

\arguments{
  \item{x}{A numeric matrix.}
  \item{y}{A numeric matrix with the same number of rows as \code{x}.}
  \item{what}{Indicates which correlations to calculate and return.  See
    value, below.}
  \item{scaled}{If TRUE, assume that the columns of \code{x} and
    \code{y} each have mean 0 and SD 1.}
  \item{corthresh}{Threshold on correlations if
    \code{what="bestpairs"}.}
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
# what = 
#  "paired":    We must have ncol(x) = ncol(y), and we then calculate 
#               the correlation for each column of x with the 
#               corresponding column of y
#
#  "bestright": We consider each column of x in turn, and find the
#               most column in y with the highest correlation,
#               returning a matrix with the max corr and the column 
#               index in y
#
#  "bestpairs": We find all pairs of columns (one in x, one in y)
#               with correlation >= corthreshold
#             
#  "all":       Calculate correlation betwen all pairs of columns
#               (one in x, one in y)
#

  If \code{what="paired"}, the return value is a vector of correlations,
  between columns of \code{x} and the corresponding column of \code{y}.
  \code{x} and \code{y} must have the same number of columns.

  If \code{what="bestright"}, we return a matrix
  of size \code{ncol(x)} by \code{2}, with the \eqn{i}th row being the maximum 
  correlation between column \eqn{i} of \code{x} and a column of
  \code{y}, and then the \code{y}-column index with that correlation.
  (In case of ties, we give the first one.)

  If \code{what="bestpairs"}, we return a matrix
  with three columns, containing all pairs of columns (with one in
  \code{x} and one in \code{y}) with correlation \eqn{\ge} \code{corthresh}.
  Each row corresponds to a column pair, and contains the correlation
  and then the \code{x}- and \code{y}-column indices.

  If \code{what="all"}, the output is a matrix of size \code{ncol(x)} by
  \code{ncol(y)}, with all correlations between columns of \code{x} and
  columns of \code{y}.
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
corbetw2mat(x, y, what="bestright")

# the same with y columns shuffled 
corbetw2mat(x, z, what="bestright")


# all pairs of columns with correlation >= 0.6
corbetw2mat(x, y, what="bestpairs", corthresh=0.6)

# the same with y columns shuffled 
corbetw2mat(x, z, what="bestpairs", corthresh=0.6)

# all correlations
corbetw2mat(x, y, what="all")
}

\seealso{ \code{\link{distee}}, \code{\link{fscale}}, \code{\link{findCommonID}} }

\keyword{array}
\keyword{univar}
\keyword{multivariate}

