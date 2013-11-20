\name{ComplexCF}
\alias{ComplexCF}
\title{
The characteristic function of the stable law
}
\description{
Theoretical characteristic function (CF) of stable law under
parametrisation \sQuote{S0} or \sQuote{S1}. See Nolan (2013) for more details.
}
\usage{
ComplexCF(t, theta, pm = 0)
}
\arguments{
  \item{t}{
    Vector of (real) numbers where the CF is evaluated ; numeric
  }
  \item{theta}{
    Vector of parameters of the stable law; vector of length 4.
  }
  \item{pm}{
    Parametrisation, an integer (0 or 1); default: \code{pm=0}(
    the Nolan \sQuote{S0} parametrisation).
  }
}
\details{
For more details about the different parametrisation of the CF, see \emph{Nolan(2013)}.
}
\value{
  Vector of complex numbers with dimension \code{length(t)}.
}
\references{
  % bibentry: nolan:2013
Nolan JP (2013).
\emph{Stable Distributions - Models for Heavy Tailed Data}.
Birkhauser, Boston.
In progress, Chapter 1 online at academic2.american.edu/\$sim\$jpnolan.
% end:bibentry:  nolan:2013
}

\seealso{
 \code{\link{jacobianComplexCF}} 
}
\examples{
# define the parameters
nt <- 10
t <- seq(0.1,3,length.out=nt)
theta <- c(1.5,0.5,1,0)
pm <- 0

# Compute the characteristic function
CF <- ComplexCF(t=t,theta=theta,pm=pm)

}

\keyword{stable-functions}


