\name{IntegrateRandomVectorsProduct}
\alias{IntegrateRandomVectorsProduct}
\title{
  Integral outer product of random vectors
}
\description{
  Computes the integral outer product of two possibly complex random
  vectors.
}
\usage{
IntegrateRandomVectorsProduct(f_fct, X, g_fct, Y, s_min, s_max,
                              subdivisions = 50,
                              method = c("Uniform", "Simpson"),
                              randomIntegrationLaw = c("norm","unif"),
                              ...)
}
\arguments{
  \item{f_fct}{
    function object with signature \code{f_fct=function(s,X)} and
    returns a matrix \eqn{ns \times nx} where nx=\code{length(X)}
    and ns=\code{length(s)}; \code{s} is the points where the integrand
    is evaluated. 
  }
  \item{X}{
    random vector where the function \code{f_fct} is evaluated. See
    Details. 
  }
  \item{g_fct}{
    function object with signature \code{g_fct=function(s,Y)} and
    returns a matrix \eqn{ns \times ny} where ny=\code{length(Y)}
    and ns=\code{length(s)}; \code{s} is the points where the integrand
    is evaluated.
  }
  \item{Y}{
    random vector where the function \code{g_fct} is evaluated. See
    Details. 
  }
  \item{s_min,s_max}{
    limits of integration. Should be finite.
  }
  \item{subdivisions}{
    maximum number of subintervals.
  }
  \item{method}{
    numerical integration rule, one of \code{"uniform"} (fast) or
    \code{"Simpson"} (more accurate quadratic rule).
  }
  \item{randomIntegrationLaw}{
    Random law pi(s) to be applied to the Random product vector, see
    Details. Choices are \code{"unif"} (uniform) and \code{"norm"}
    (normal distribution).
  }
  \item{\dots}{
    other arguments to pass to random integration law. Mainly, the mean
    (\code{mu}) and standard deviation (\code{sd}) of the normal law.
  }
}
\details{

  The function computes the \eqn{nx \times ny} matrix \eqn{C =
  \int_{s_{min}}^{s_{max}} f_s(X) g_s(Y) \pi(s) ds}, such as the one
  used in the objective function of the Cgmm method. This is essentially
  an outer product with with multiplication replaced by integration.

  There is no function in R to compute vectorial integration and
  computing \eqn{C} element by element using \code{integrate} may
  be very slow when \code{length(X)} (or \code{length(y)}) is large.

  The function allows complex vectors as its integrands.
}
\value{
  an \eqn{nx \times ny} matrix \eqn{C} with elements:
  \deqn{c_{ij} =
    \int_{s_{min}}^{s_{max}} f_s(X_i) g_s(Y_j) \pi(s) ds
    .
  }
}
\examples{
## Define the integrand
f_fct <- function(s, x) {
    sapply(X = x, FUN = sampleComplexCFMoment, t = s, theta = theta)
}
f_bar_fct <- function(s, x) Conj(f_fct(s, x))

## Function specific arguments
theta <- c(1.5, 0.5, 1, 0)
set.seed(345)
X <- rstable(3, 1.5, 0.5, 1, 0)
s_min <- 0;
s_max <- 2
numberIntegrationPoints <- 10
randomIntegrationLaw <- "norm"

Estim_Simpson <-
    IntegrateRandomVectorsProduct(f_fct, X, f_bar_fct, X, s_min, s_max, 
                                  numberIntegrationPoints, 
                                  "Simpson", randomIntegrationLaw)
Estim_Simpson
}
\keyword{general-functions}
