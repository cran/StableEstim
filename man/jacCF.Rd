\name{CF Jacobian}
\alias{jacobianComplexCF}
\title{
Jacobian of the characteristic function of the stable law.
}
\description{
  Numeric jacobian of the  characteristic function (CF) as a function of
  the parameter \eqn{\theta} evaluated at a specific point(vector)
  \code{t} and a given value \eqn{\theta}.
}
\usage{
jacobianComplexCF(t, theta, pm = 0)
}
\arguments{
  \item{t}{
    Vector of (real) numbers where the jacobian of the CF is evaluated ;numeric
  }
  \item{theta}{
    Given value of the parameters where the jacobian of the CF is evaluated; vector of length 4.
  }
  \item{pm}{
    Parametrisation, an integer (0 or 1); default: \code{pm=0} (
    the Nolan \sQuote{S0} parametrisation).
  }
}
\details{
  The numerical derivation is obtained by a call to the function
  \code{jacobian} from the package \pkg{numDeriv}. We
  have set up its arguments by default and the user is not allowerd
  to modify them.
}
\value{ 
  Returns a matrix \code{length(t)} \eqn{\times } 4 of complex numbers.
}

\seealso{
 \code{\link{ComplexCF}} 
}
\examples{
# define the parameters
nt <- 10
t <- seq(0.1,3,length.out=nt)
theta <- c(1.5,0.5,1,0)
pm <- 0

# Compute the jacobian of the characteristic function
jack_CF <- jacobianComplexCF(t=t,theta=theta,pm=pm)
}

\keyword{stable-functions}
