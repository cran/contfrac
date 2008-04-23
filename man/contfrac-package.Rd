\name{contfrac-package}
\alias{contfrac-package}
\alias{contfrac}
\docType{package}
\title{
Continued fractions
}
\description{
Various utilities for manipulating continued fractions
}
\details{
\tabular{ll}{
Package: \tab contfrac\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-04-04\cr
License: \tab GPL \cr
}

}
\author{
Robin K. S. Hankin

Maintainer: <r.hankin@noc.soton.ac.uk>
}
\references{
  \item W. H. Press, B. P. Flannery, S. A. Teukolsky, and
  W. T. Vetterling 1992.  \emph{Numerical recipes 3rd edition: the art of
    scientific computing}.  Cambridge University Press; section 5.2
  \dQuote{Evaluation of continued fractions}
}
\keyword{ package }
\examples{
  # Some convergents of pi:
  jj <- convergents(c(3,7,15,1,292))
  jj$A / jj$B - pi


  # An identity of Euler's:
  jj <- GCF(a=seq(from=2,by=2,len=30), b=seq(from=3,by=2,len=30), b0=1) 
  jj - 1/(exp(0.5)-1)   # should be small


  # Now a continued fraction representation of tan(z):
  tan_cf <- function(z,n=14){ GCF(c(z,rep(-z^2,n-1)), seq(from=1,by=2,len=n)) }

  tan_cf(1+1i) - tan(1+1i)   # should be small
}
