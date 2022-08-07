\name{Estim}
\alias{Estim}
\title{Estimate parameters of stable laws}
\description{

  Estimates the four parameters of stable distributions using one of the
  methods implemented in \pkg{StableEstim}. This is the main user-level
  function but the individul methods are available also as separate
  functions.  

}
\usage{
Estim(EstimMethod = c("ML", "GMM", "Cgmm","Kout"), data, theta0 = NULL,
      ComputeCov = FALSE, HandleError = TRUE, ...)
}
\arguments{
  \item{EstimMethod}{
    Estimation method to be used, one of \code{"ML"} (maximum
    likelihood, default), \code{"GMM"} (generalised method of moment with finite
    moment conditions), \code{"Cgmm"} (GMM with continuum moment
    conditions), and \code{"Kout"} (Koutrouvelis regression method).
  }
  \item{data}{
    Data used to perform the estimation, a numeric vector.
  }
  \item{theta0}{
    Initial values for the 4 parameters. If \code{NULL} (default),
    initial values are computed using the fast Kogon-McCulloch method,
    see \code{\link{IGParametersEstim}}; vector of length 4.
  }
  \item{ComputeCov}{
    Logical flag: if \code{TRUE}, the asymptotic covariance matrix (4x4)
    is computed (except for the Koutrouvelis method).
  }
  \item{HandleError}{
    Logical flag: if \code{TRUE} and if an error occurs during the
    estimation procedure, the computation will carry on and \code{NA}
    will be returned. Useful for Monte Carlo simulations, see
    \code{\link{Estim_Simulation}}.
  }
  \item{\dots}{
    Other arguments to be passed to the estimation function, such as the
    asymptotic confidence level, see Details.
  }
}
\details{

  \code{Estim} is the main estimation function in package
  \pkg{StableEstim}. 

  This function should be used in priority for estimation purpose as it
  provides more information about the estimator. However, user needs to
  pass the appropriate parameters to the selected method in
  \code{\dots}. See the documentation of the selected method.

  \bold{Asymptotic Confidence Intervals}:
  The \emph{normal} asymptotic confidence intervals (CI) are computed.
  The user can set the \emph{level} of confidence by inputing the
  \code{level} argument (in the \code{"\dots"}); default
  \code{level=0.95}. The theoretical justification for asymptotic normal
  CI can be found in the references for the individual methods. Note the
  CI's are not computed for the Koutrouvelis regression method.

}
\value{
  an object of class \code{Estim}, see \code{\link{Estim-class}} for
  more details
}

\seealso{
  
  \code{\link{CgmmParametersEstim}},
  \code{\link{GMMParametersEstim}},
  \code{\link{MLParametersEstim}},
  \code{\link{KoutParametersEstim}} for the individual estimation
  methods;

  \code{\link{IGParametersEstim}} for fast computation of initial
  values. 
  
}
% CHECK
\examples{
## general inputs
theta <- c(1.45, 0.55, 1, 0)
pm <- 0
set.seed(2345)
x <- rstable(200, theta[1], theta[2], theta[3], theta[4], pm)

objKout <- Estim(EstimMethod = "Kout", data = x, pm = pm, 
                     ComputeCov = FALSE, HandleError = FALSE, 
                     spacing = "Kout")
}
\keyword{Estim-functions}
