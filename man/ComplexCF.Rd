\name{ComplexCF}
\alias{ComplexCF}
\title{
  Compute the characteristic function of stable laws
}
\description{
  Theoretical characteristic function (CF) of stable laws under
  parametrisation \sQuote{S0} or \sQuote{S1}. See Nolan (2013) for more
  details.
}
\usage{
ComplexCF(t, theta, pm = 0)
}
\arguments{
  \item{t}{
    vector of (real) numbers where the CF is evaluated; numeric.
  }
  \item{theta}{
    vector of parameters of the stable law; vector of length 4.
  }
  \item{pm}{
    parametrisation, an integer (0 or 1); default: \code{pm = 0}
    (Nolan's \sQuote{S0} parametrisation).
  }
}
\details{
  For more details about the different parametrisation of the CF, see
  \emph{Nolan(2012)}.
}
\value{
  vector of complex numbers with dimension \code{length(t)}.
}
\references{
  \insertRef{nolan:2012}{StableEstim}
}

\seealso{
  \code{\link{jacobianComplexCF}} 
}
\examples{
## define the parameters
nt <- 10
t <- seq(0.1, 3, length.out = nt)
theta <- c(1.5, 0.5, 1, 0)
pm <- 0

## Compute the characteristic function
CF <- ComplexCF(t = t, theta = theta, pm = pm)
CF
}
\keyword{stable-functions}
