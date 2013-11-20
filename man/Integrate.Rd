\name{IntegrateRandomVectorsProduct}
\alias{IntegrateRandomVectorsProduct}
\title{
 Integrate Ranndom vectors product
}
\description{
Computes the integral of two random vector.
}
\usage{
IntegrateRandomVectorsProduct(f_fct, X, g_fct, Y, s_min, s_max,
                              subdivisions = 50,
                              method = c("Uniform", "Simpson"),
                              randomIntegrationLaw = c("norm","unif"),
                              ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f_fct}{
    Function object with signature \code{f_fct=function(s,X)} and
    returns a matrix \eqn{ns \times nx} where nx=\code{length(X)}
    and ns=\code{length(s)}; \code{s} is the points where the integrand
    is evaluated. 
  }
  \item{X}{
    Random vector where the function \code{f_fct} is evaluated. See details
  }
  \item{g_fct}{
    Function object with signature \code{g_fct=function(s,Y)} and
    returns a matrix \eqn{ns \times ny} where ny=\code{length(Y)}
    and ns=\code{length(s)}; \code{s} is the points where the integrand
    is evaluated.
  }
  \item{Y}{
    Random vector where the function \code{g_fct} is evaluated. See details
  }
  \item{s_min,s_max}{
    Limits of integration. Should be finite
  }
  \item{subdivisions}{
    Maximum number of subintervals.
  }
  \item{method}{
    Numerical integration rule. User ca choose between
    \code{"uniform"} (fast) or \code{"Simpson"} (more accurate quadratic rule). 
  }
  \item{randomIntegrationLaw}{
    Random law pi(s) to be applied to the Random product vector, see
    details. Choices are \code{"unif"} (uniform) and \code{"norm"}
    (normal distribution).
  }
  \item{\dots}{
    Other arguments to pass to random integration law. Mainly, the mean
    (\code{mu}) and standard deviation (\code{sd}) of the normal law.
  }
}
\details{
  %\cr\cr
  The function computes the \eqn{nx \times ny} matrix
  \eqn{C=\int_{s_{min}}^{s_{max}} f_s(X) g_s(Y) \pi(s) ds} such as the one
  used in the objective function of the Cgmm method.
  %\cr\cr
  There is no function in R to compute vectorial integration and
  computing \eqn{C} element by element using \code{integrate} may
  be very slow when \code{length(X)} (or \code{length(y)}) is large.
  %\cr\cr
  The function allows complex vectors as its integrands.
}
\value{
  Returns a \eqn{nx \times ny} matrix \eqn{C} with elements:
  \deqn{c_{ij}=\int_{s_{min}}^{s_{max}} f_s(X_i) g_s(Y_j) \pi(s) ds} 
}

\examples{
## Define the integrand
f_fct <- function(s,x){sapply(X=x,
                              FUN=sampleComplexCFMoment,
                              t=s,theta=theta)
                   }
f_bar_fct <- function(s,x){Conj(f_fct(s,x))}

## Function specific arguments
theta <- c(1.5,0.5,1,0)
set.seed(345);X=rstable(3,1.5,0.5,1,0)
s_min=0;s_max=2
numberIntegrationPoints=10
randomIntegrationLaw="norm"

Estim_Simpson<- IntegrateRandomVectorsProduct(f_fct,X,f_bar_fct,X,s_min,s_max,
                                              numberIntegrationPoints,
                                              "Simpson",randomIntegrationLaw)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{general-functions}% __ONLY ONE__ keyword per line
