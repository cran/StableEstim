\name{StableEstim-package}
\alias{StableEstim-package}
\docType{package}
\title{
  Stable law estimation functions
}
\description{
  A collection of methods to estimate the four parameters of stable
  laws. The package also provides functions to compute the
  characteristic function and tools to run Monte Carlo simulations.
}
\details{
  The main functions of the package are briefly described below:
  \describe{
    \item{main function:}{\code{\link{Estim}} is the most useful
      function of the package. It estimates of the parameters
      and the asymptotic properties of the estimators.}

    \item{estimation function:}{
      the methods provided so far are the maximum-likelihood
      (\code{\link{MLParametersEstim}}), the generalised method of
      moment with finite (\code{\link{GMMParametersEstim}}) or
      continuum (\code{\link{CgmmParametersEstim}}) moment conditions,
      the iterative Koutrouvelis regression method
      (\code{\link{KoutParametersEstim}}) and the fast Kogon-McCulloch
      method used for first guess estimation
      (\code{\link{IGParametersEstim}}).
    }
    \item{characteristic function:}{the characteristic function
      (\code{\link{ComplexCF}})
      and its Jacobian (\code{\link{jacobianComplexCF}})
      can be computed and will return a vector (respectively a matrix)
      of complex numbers.}
    \item{Monte Carlo simulation}{
      \code{\link{Estim_Simulation}} is a tool to run Monte Carlo
      simulations with flexible options to select the estimation method,
      the Monte Carlo control parameters, compute statistical summaries
      or save results to a file.
    }
  }
}
\author{
Tarak Kharrat, Georgi N. Boshnakov
}
\note{
  Version 1 of this package had a somewhat restricted license since
  it needed package \pkg{akima} in some computations.

  In version 2 of the package we implemented a 2D interpolation routine
  and removed the dependency on \pkg{akima}.  Therefore,
  \pkg{StableEstim} is now under GPL license.  The package is related to
  upcoming work by the authors where the different methods are compared
  using MC simulations.

}
\seealso{
  \code{fBasics:::.mleStableFit},
  \code{fBasics:::.qStableFit}
  
  package \pkg{stabledist}
}
\references{% bibentry:all
Carrasco M and Florens J (2000).
``Generalization of GMM to a continuum of moment conditions.''
\emph{Econometric Theory}, \bold{16}(06), pp. 797--834.

Carrasco M and Florens J (2002).
``Efficient GMM estimation using the empirical characteristic function.''
\emph{IDEI Working Paper}, \bold{140}.

Carrasco M and Florens J (2003).
``On the asymptotic efficiency of GMM.''
\emph{IDEI Working Paper}, \bold{173}.

Carrasco M, Chernov M, Florens J and Ghysels E (2007).
``Efficient estimation of general dynamic models with a continuum of moment conditions.''
\emph{Journal of Econometrics}, \bold{140}(2), pp. 529--573.

Carrasco M, Florens J and Renault E (2007).
``Linear inverse problems in structural econometrics estimation based on spectral decomposition and regularization.''
\emph{Handbook of econometrics}, \bold{6}, pp. 5633--5751.

Carrasco M and Kotchoni R (2010).
``Efficient estimation using the characteristic function.''
Mimeo. University of Montreal.

Nolan J (2001).
``Maximum likelihood estimation and diagnostics for stable distributions.''
\emph{L'evy processes: theory and applications}, pp. 379--400.


\insertRef{nolan:2012}{StableEstim}


Hansen LP (1982).
``Large sample properties of generalized method of moments estimators.''
\emph{Econometrica: Journal of the Econometric Society}, pp. 1029--1054.

Hansen LP, Heaton J and Yaron A (1996).
``Finite-sample properties of some alternative GMM estimators.''
\emph{Journal of Business & Economic Statistics}, \bold{14}(3), pp. 262--280.

Feuerverger A and McDunnough P (1981).
``On efficient inference in symmetric stable laws and processes.''
\emph{Statistics and Related Topics}, \bold{99}, pp. 109--112.

Feuerverger A and McDunnough P (1981).
``On some Fourier methods for inference.''
\emph{Journal of the American Statistical Association}, \bold{76}(374), pp. 379--387.

Schmidt P (1982).
``An improved version of the Quandt-Ramsey MGF estimator for mixtures of normal distributions and switching regressions.''
\emph{Econometrica: Journal of the Econometric Society}, pp. 501--516.

Besbeas P and Morgan B (2008).
``Improved estimation of the stable laws.''
\emph{Statistics and Computing}, \bold{18}(2), pp. 219--231.
% end:bibentry:all
}
\keyword{package}
