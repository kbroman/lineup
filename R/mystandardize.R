mystandardize <-
function(x)
{
  if(!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  matrix(.C("R_mystandardize",
            as.integer(n),
            as.integer(p),
            x=as.double(x),
            PACKAGE="lineup",
            NAOK=TRUE)$x, ncol=p)
}
