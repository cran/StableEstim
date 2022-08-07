\name{CgmmParametersEstim}
\alias{CgmmParametersEstim}
\title{
  Estimate parameters of stable laws using a Cgmm method
}
\description{
  Estimate the four parameters of stable laws using generalised method
  of moments based on a continuum of complex moment conditions (Cgmm)
  due to Carrasco and Florens.  Those moments are computed by matching
  the characteristic function with its sample counterpart. The resulting
  (ill-posed) estimation problem is solved by a regularisation
  technique.
}
\usage{
CgmmParametersEstim(x, type = c("2S", "IT", "Cue"), alphaReg = 0.01,
                    subdivisions = 50,
                    IntegrationMethod = c("Uniform", "Simpson"),
                    randomIntegrationLaw = c("unif", "norm"),
                    s_min = 0, s_max = 1,
                    theta0 = NULL,
                    IterationControl = list(),
                    pm = 0, PrintTime = FALSE,...)
}
\arguments{
  \item{x}{
    Data used to perform the estimation: a vector of length n.
  }
  \item{type}{
    Cgmm algorithm: \code{"2S"} is the two steps GMM proposed by
    Hansen(1982).  \code{"Cue"} and \code{"IT"} are respectively the
    continuous updated and the iterative GMM proposed by Hansen, Eaton
    et Yaron (1996) and adapted to the continuum case.
  }
  \item{alphaReg}{
    Value of the regularisation parameter; numeric, default = 0.01.
  }
  \item{subdivisions}{
    Number of subdivisions used to compute the different
    integrals involved in the computation of the objective function (to
    minimise); numeric.
  }
  \item{IntegrationMethod}{
    Numerical integration method to be used to approximate the
    (vectorial) integrals. Users can choose between \code{"Uniform"}
    discretization or the \code{"Simpson"}'s rule (the 3-point
    Newton-Cotes quadrature rule).
  }
  \item{randomIntegrationLaw}{
    Probability measure associated to the Hilbert space spanned by the
    moment conditions. See Carrasco and Florens (2003) for more details.
  }
  \item{s_min,s_max}{
    Lower and Upper bounds of the interval where the
    moment conditions are considered; numeric.
  }
  \item{theta0}{
    Initial guess for the 4 parameters values: vector of length 4.
  }
  \item{IterationControl}{
    Only used with \code{type = "IT"} or \code{type = "Cue"} to control the
    iterations, see Details.
  }
  \item{pm}{
    Parametrisation, an integer (0 or 1); default: \code{pm = 0} (Nolan's
    \sQuote{S0} parametrisation).
  }
  \item{PrintTime}{
    Logical flag; if set to TRUE, the estimation duration is printed out to
    the screen in a readable format (h/min/sec).
  }
  \item{\dots}{
    Other arguments to be passed to the optimisation function and/or to
    the integration function.
  }
}
\details{
  \bold{The moment conditions}
  %\cr\cr
  The moment conditions are given by:
  \deqn{g_t(X,\theta)=g(t,X;\theta)= e^{itX} - \phi_{\theta}(t)}
  If one has a sample \eqn{x_1,\dots,x_n} of i.i.d realisations of the
  same random variable \eqn{X}, then: 
  \deqn{\hat{g}_n(t,\theta)  = \frac{1}{n}\sum_{i=1}^n g(t,x_i;\theta) =  \phi_n(t) -\phi_\theta(t),}
  where \eqn{\phi_n(t)} is the eCF associated with the sample
  \eqn{x_1,\dots,x_n}, defined by \eqn{\phi_n(t)= \frac{1}{n}
    \sum_{j=1}^n e^{itX_j}}.
  %\cr\cr
  \bold{Objective function}
  %\cr\cr
  
  Following \insertCite{@@carrasco2007efficientCont, Proposition 3.4;textual}{StableEstim},
  the objective function to minimise is given by:
  \deqn{obj(\theta)=\overline{\underline{v}^{\prime}}(\theta)[\alpha_{Reg} \mathcal{I}_n+C^2]^{-1}\underline{v}(\theta)}
  where:
  \describe{
    \item{
      \eqn{\underline{v} = [v_1,\ldots,v_n]^{\prime}};}{\eqn{v_i(\theta)
	= \int_I \overline{g_i}(t;\hat{\theta}^1_n) \hat{g}(t;\theta) \pi(t) dt}.
    }
    \item{\eqn{I_n}}{is the identity matrix of size \eqn{n}.}
    \item{\eqn{C}}{is a \eqn{n \times n} matrix with \eqn{(i,j)}th
      element given by
      \eqn{c_{ij} = \frac{1}{n-4}\int_I
	\overline{g_i}(t;\hat{\theta}^1_n) g_j(t;\hat{\theta}^1_n)
	\pi(t) dt}.
    }
  }
  To compute \eqn{C} and \eqn{v_i()} we will use the function
  \code{\link{IntegrateRandomVectorsProduct}}.
  %\cr\cr
  \bold{The IterationControl}
  %\cr\cr
  If \code{type = "IT"} or \code{type = "Cue"}, the user can control
  each iteration using argument \code{IterationControl}, which should be
  a \code{list} which contains the following elements:
  \describe{
    \item{\code{NbIter}:}{maximum number of iterations.
      The loop stops when \code{NBIter} is reached; default = 10.
    }
    \item{\code{PrintIterlogical}:}{if set to TRUE the values of the
      current parameter estimates are printed to the screen at each
      iteration; default = TRUE.
    }
    \item{\code{RelativeErrMax}:}{the loop stops if the relative error
      between two consecutive estimation steps is smaller then
      \code{RelativeErrMax}; default = 1e-3.
    }
  }
}
\value{
  a list with the following elements:
  \item{Estim}{output of the optimisation function,}
  \item{duration}{estimation duration in numerical format,}
  \item{method}{\code{character} describing the method used.}
}

\references{
  
  \insertRef{carrasco2000generalization}{StableEstim}

  \insertRef{carrasco2002efficient}{StableEstim}

  \insertRef{carrasco2003asymptotic}{StableEstim}

  \insertRef{carrasco2007efficientCont}{StableEstim}

  \insertRef{carrasco2010efficient}{StableEstim}
}
\note{
  \code{nlminb} as used to minimise the Cgmm objective function.
}
\seealso{
  \code{\link{Estim}},
  \code{\link{GMMParametersEstim}},
  \code{\link{IntegrateRandomVectorsProduct}}
}
\examples{
## general inputs
theta <- c(1.45, 0.55, 1, 0)
pm <- 0
set.seed(2345)
x <- rstable(50, theta[1], theta[2], theta[3], theta[4], pm)

## GMM specific params
alphaReg <- 0.01
subdivisions <- 20
randomIntegrationLaw <- "unif"
IntegrationMethod <- "Uniform"

## Estimation
twoS <- CgmmParametersEstim(x = x, type = "2S", alphaReg = alphaReg, 
                          subdivisions = subdivisions, 
                          IntegrationMethod = IntegrationMethod, 
                          randomIntegrationLaw = randomIntegrationLaw, 
                          s_min = 0, s_max = 1, theta0 = NULL, 
                          pm = pm, PrintTime = TRUE)
twoS
}
\keyword{Estim-functions}
