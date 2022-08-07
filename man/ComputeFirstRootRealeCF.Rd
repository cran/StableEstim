\name{ComputeFirstRootRealeCF}
\alias{ComputeFirstRootRealeCF}
\title{
  First root of the empirical characteristic function
}
\description{
  Computes the first root of the real part of the empirical
  characteristic function. 
}
\usage{
ComputeFirstRootRealeCF(x, ..., tol = 0.001, maxIter = 100,
                        lowerBand = 1e-04, upperBand = 30)
}
\arguments{
  \item{x}{
    data used to perform the estimation: vector of length n.
  }
  \item{\dots}{
    other arguments to pass to the optimisation function.
  }
  \item{tol}{
    tolerance to accept the solution; default = 1e-3.
  }
  \item{maxIter}{
    maximum number of iteration in the Welsh algorithm; default = 100.
  }
  \item{lowerBand}{
    lower band of the domain where the graphical seach is performed;
    default = 1e-4.
  }
  \item{upperBand}{
    Lower band of the domain where the graphical seach is performed;
    default = 30.
  }
}
\details{
  The Welsh algorithm is first applied. If it fails to provide a
  satisfactory value (\code{< tol}), a graphical/ numerical approach is
  used. We first plot the real part of the eCF vs t in order to
  determine the first zero directly and use it as the initial guess of a
  numerical minimisation routine.
}
\value{
  \code{numeric}: first zero of the real part of the eCF.
}
\references{
  \insertRef{welsh1986implementing}{StableEstim}
}
\seealso{
  \code{\link{ComplexCF}}
}
\examples{
set.seed(345)
x <- rstable(500, 1.5, 0.5)
ComputeFirstRootRealeCF(x)
}
\keyword{stable-functions}
