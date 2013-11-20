\name{Estim}
\alias{Estim}
\title{Main estimation function}
\description{
Main estimation function which computes all the information about the
estimator (and its asymptotic properties). It allows the user to choose
the preferred method and several related options.
}
\usage{
Estim(EstimMethod = c("ML", "GMM", "Cgmm","Kout"), data, theta0 = NULL,
      ComputeCov = FALSE, HandleError = TRUE, ...)
}
\arguments{
  \item{EstimMethod}{
    Estimation method to be used. User can choose between ML,
    GMM, Cgmm or Koutrouvelis regression method.
  }
  \item{data}{
    Data used to perform the estimation: vector of length n.
  }
  \item{theta0}{
    Initial guess for the 4 parameters values: If \code{NULL} the
    Kogon-McCulloch method is called, see
    \code{\link{IGParametersEstim}}; vector of length 4.
  }
  \item{ComputeCov}{
    Logical flag: If set to TRUE, the asymptotic covariance matrix (4x4) is
    computed (except for the Koutrouvelis method).
  }
  \item{HandleError}{
    Logical flag: If set to TRUE and if an error occurs during the estimation
    procedure, the computation will carry on and NA will be returned. Useful
    for Monte Carlo simulations, see \code{\link{Estim_Simulation}}
  }
  \item{\dots}{
    Other arguments to be passed to the estimation function or the
    asymptotic confidence level. See details.
  }
}
\details{
  %\cr\cr
  This function should be used in priority for estimation purpose as it
  provides more information about the estimator. However, user needs to
  pass the appropriate parameters to the selected method in
  \code{\dots}. See the documentation of the selected method.
  %\cr\cr
  \bold{Asymptotic Confidence Interval}:
  %\cr\cr
  The \emph{normal} asymptotic confidence interval (C.I) are computed.
  The user can set the \emph{level} of confidence by inputting the
  \code{level} argument (in the \code{\dots}); default
  \code{level}=0.95. The theoretical justification for asymptotic
  normal C.I could be find in the references of each method. Note the
  C.I are not computed for the Koutrouvelis regression method.
}
\value{
  Returns an object of class \code{Estim}. See \code{\link{Estim}} for
  more details. 
}

\seealso{
 \code{\link{CgmmParametersEstim}}, \code{\link{GMMParametersEstim}},
 \code{\link{MLParametersEstim}}, \code{\link{KoutParametersEstim}}
}

% CHECK
\examples{
## general inputs
theta <- c(1.45,0.55,1,0)
pm <- 0
set.seed(2345)
x <- rstable(200,theta[1],theta[2],theta[3],theta[4],pm)

objKout <- Estim(EstimMethod="Kout",data=x,pm=pm,
                     ComputeCov=FALSE,HandleError=FALSE,
                     spacing="Kout")

}

\keyword{Estim-functions}

