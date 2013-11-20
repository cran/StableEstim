\name{CgmmParametersEstim}
\alias{CgmmParametersEstim}
\title{
The generalised method of moment with a continuum of moment conditions
}
\description{
Cgmm method of Carrasco and Florens to estimate the 4 parameters of
stable law based on a continuum of complex moment conditions (the modulus).
Those moments are computed by matching the characteristic
function with its sample counterpart. The resulting (ill-posed) estimation
problem is solved by a regularisation technique.
}
\usage{
CgmmParametersEstim(x, type = c("2S", "IT", "Cue"), alphaReg = 0.01,
                    subdivisions=50,
                    IntegrationMethod = c("Uniform", "Simpson"),
                    randomIntegrationLaw = c("unif", "norm"),
                    s_min=0,s_max = 1,
                    theta0 = NULL,
                    IterationControl = list(),
                    pm = 0, PrintTime = FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Data used to perform the estimation: a vector of length n.
  }
  \item{type}{
    Cgmm algorithm: "2S" is the two steps GMM proposed by Hansen(1982)
    and the "Cue" and "IT" are respectively the continuous updated and the
    iterative GMM proposed by Hansen, Eaton et Yaron (1996) and adapted
    to the continuum case.
  }
  \item{alphaReg}{
   Value of the regularisation parameter; numeric, default=0.01
 }
 \item{subdivisions}{
   Number of subdivision used to compute the different
   integrals involved in the computation of the objective function (to
   minimise); numeric
 }
 \item{IntegrationMethod}{
   Numerical integration method to be used to approximate the (vectorial)
   integrals. User can choose between the "Uniform" discretization or the
   "Simpson"'s rule (the 3-point Newton-Cotes quadrature rule).
 }
 \item{randomIntegrationLaw}{
   Probability measure associated to the Hilbert space spanned by the
   moment conditions. See Carrasco and Florens (2003) for more details.
 }
 \item{s_min,s_max}{Lower and Upper bounds of the interval where the
   moment conditions are considered; numeric
 }
 \item{theta0}{
   Initial guess for the 4 parameters values: vector of length 4
 }
 \item{IterationControl}{
   Only used with \code{type="IT"} or \code{type="Cue"} to control the
   iterations. See details
 }
 \item{pm}{
   Parametrisation, an integer (0 or 1); default: \code{pm=0} (the Nolan
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
  If one has a sample \eqn{x_1,\dots,x_n} of i.i.d realisations of the same random variable \eqn{X}, then:
  \deqn{\hat{g}_n(t,\theta)  = \frac{1}{n}\sum_{i=1}^n g(t,x_i;\theta) =  \phi_n(t) -\phi_\theta(t)}
  where \eqn{\phi_n(t)} is the eCF associated to the sample
  \eqn{x_1,\dots,x_n} and defined by \eqn{\phi_n(t)= \frac{1}{n}
    \sum_{j=1}^n e^{itX_j}}
  %\cr\cr
  \bold{Objective function}
  %\cr\cr
  Following Carrasco and al (Proposition 3.4, 2007), the objective
  function to minimise is given by:
  \deqn{obj(\theta)=\overline{\underline{v}^{\prime}}(\theta)[\alpha_{Reg} \mathcal{I}_n+C^2]^{-1}\underline{v}(\theta)}
 where:
  \describe{
    \item{\eqn{\underline{v}=[v_1,\ldots,v_n]^{\prime}} ;}{\eqn{v_i(\theta)=\int_I \overline{g_i}(t;\hat{\theta}^1_n) \hat{g}(t;\theta) \pi(t) dt}}
    \item{\eqn{I_n}}{is the identity matrix of size \eqn{n}}
    \item{\eqn{C}}{is a \eqn{n \times n} matrix with (i,j) elements given by
      \eqn{c_{ij}=\frac{1}{n-4}\int_I \overline{g_i}(t;\hat{\theta}^1_n) g_j(t;\hat{\theta}^1_n) \pi(t) dt}}
    }
    To compute C and \eqn{v_i()} we will use the function
    \code{\link{IntegrateRandomVectorsProduct}}.
    %\cr\cr
    \bold{The IterationControl}
    %\cr\cr
    If \code{type="IT"} or \code{type="Cue"} is selected, user can control each
    iteration by setting up the \code{list} \code{IterationControl} which
    contains the following elements:
    \describe{
      \item{\code{NbIter}:}{maximum number of iteration. The loop stops when \code{NBIter}
	is reached; default=10}
      \item{\code{PrintIterlogical}:}{if set to TRUE the value of the current
	parameter estimation is printed to the screen at each iteration; default=TRUE}
      \item{\code{RelativeErrMax}:}{the loop stops if the relative error between two
	consecutive estimation is smaller then
	\code{RelativeErrMax}; default=1e-3}
    }
  }
  \value{
    Returns a list with the following elements:
    \item{Estim}{output of the optimisation function}
    \item{duration}{estimation duration in numerical format}
    \item{method}{\code{character} describing the method used}
  }

\references{
  % bibentry: carrasco2000generalization
Carrasco M and Florens J (2000).
``Generalization of GMM to a continuum of moment conditions.''
\emph{Econometric Theory}, \bold{16}(06), pp. 797--834.
% end:bibentry:  carrasco2000generalization

  % bibentry: carrasco2002efficient
Carrasco M and Florens J (2002).
``Efficient GMM estimation using the empirical characteristic function.''
\emph{IDEI Working Paper}, \bold{140}.
% end:bibentry:  carrasco2002efficient

  % bibentry: carrasco2003asymptotic
Carrasco M and Florens J (2003).
``On the asymptotic efficiency of GMM.''
\emph{IDEI Working Paper}, \bold{173}.
% end:bibentry:  carrasco2003asymptotic

  % bibentry: carrasco2007efficientCont
Carrasco M, Chernov M, Florens J and Ghysels E (2007).
``Efficient estimation of general dynamic models with a continuum of moment conditions.''
\emph{Journal of Econometrics}, \bold{140}(2), pp. 529--573.
% end:bibentry:  carrasco2007efficientCont

  % bibentry: carrasco2010efficient
Carrasco M and Kotchoni R (2010).
``Efficient estimation using the characteristic function.''
Mimeo. University of Montreal.
% end:bibentry:  carrasco2010efficient
}

\note{
  \code{nlminb} as used to minimise the Cgmm objective function.
}

\seealso{
\code{\link{Estim}}, \code{\link{IntegrateRandomVectorsProduct}}
}
\examples{ %CHECK
## general inputs
theta <- c(1.45,0.55,1,0)
pm <- 0
set.seed(2345);x <- rstable(50,theta[1],theta[2],theta[3],theta[4],pm)

## GMM specific params
alphaReg=0.01
subdivisions=20
randomIntegrationLaw="unif"
IntegrationMethod="Uniform"

## Estimation
twoS <- CgmmParametersEstim(x=x,type="2S",alphaReg=alphaReg,
                          subdivisions=subdivisions,
                          IntegrationMethod=IntegrationMethod,
                          randomIntegrationLaw=randomIntegrationLaw,
                          s_min=0,s_max=1,theta0=NULL,
                          pm=pm,PrintTime=TRUE)
}
\keyword{Estim-functions}
