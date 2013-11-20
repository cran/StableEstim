\name{expect_almost_equal}
\alias{expect_almost_equal}
\title{
  test approximate equality
}
\description{
  tests the approximate equality of 2 objects. Useful for running tests.
}
\usage{
expect_almost_equal(x, y, tolExpect = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    First object.
}
  \item{y}{
    Second object
}
  \item{tolExpect}{
    Tolerance accepted; default=1e-3
  }
}
\details{
  This function works with the \code{expect_that} function from package
  \code{testhat} to test equality between 2 objects with a given
  tolerance. It used particularly for testing functions output. See the
  CF examples in the Examples folder.
}

\seealso{
  \code{expect_that},\code{testthat}
}
\examples{
x=1.1 ; y=1.5
expect_almost_equal(x,y,1) ## pass
#expect_almost_equal(x,y,0.3) -> fails
}

\keyword{general-functions}
