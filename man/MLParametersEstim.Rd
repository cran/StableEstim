\name{MLParametersEstim}
\alias{MLParametersEstim}
\title{
  Maximum likelihood (ML) method
}
\description{
  Uses the numerical ML approach described by Nolan to estimate the 4
  parameters of stable law. The method may be slow for large sample size
  due to the use of numerical optimisation routine.
}
\usage{
MLParametersEstim(x, theta0 = NULL, pm = 0, PrintTime = FALSE, ...)
}
\arguments{
  \item{x}{
    data used to perform the estimation: vector of length n.
  }
  \item{theta0}{
    initial guess for the 4 parameters values: If \code{NULL}, the
    Kogon-McCulloch method is called, see \code{\link{IGParametersEstim}};
    a vector of length 4.
  }
  \item{pm}{
    parametrisation, an integer (0 or 1); default: \code{pm=0}
    (Nolan's \sQuote{S0} parametrisation).
  }
  \item{PrintTime}{
    logical flag; if set to TRUE, the estimation duration is printed out
    to the screen in a readable format (h/min/sec).
  }
  \item{\dots}{
    Other argument to be passed to the optimisation function.
  }
}
\details{
  The function performs the minimisation of the numerical (-)log-density
  of stable laws computed by function \code{dstable} from package
  \pkg{stabledist}.

  After testing several optimisation routines, we have found out that
  the \code{"L-BFGS-B"} algorithm performs better with the ML method
  (faster, more accurate).
}
\value{
  a list with the following elements:
  \item{Estim }{output of the optimisation function,}
  \item{duration}{estimation duration in a numerical format,}
  \item{method}{\code{character} describing the method used.}
}
\references{
  % bibentry: nolan2001maximum
Nolan J (2001).
``Maximum likelihood estimation and diagnostics for stable distributions.''
\emph{L'evy processes: theory and applications}, pp. 379--400.
% end:bibentry:  nolan2001maximum
}
\seealso{
  \code{\link{Estim}}
}
%CHECK
\examples{
theta <- c(1.5, 0.4, 1, 0)
pm <- 0
## 50 points does not give accurate estimation
## but it makes estimation fast for installation purposes
## use at least 200 points to get decent results.
set.seed(1333)
x <- rstable(50, theta[1], theta[2], theta[3], theta[4], pm)

## This example takes > 30 sec hence commented out
\dontrun{
  ML <- MLParametersEstim(x = x, pm = pm, PrintTime = TRUE)
}
## see the Examples folder for more examples.
}
\keyword{Estim-functions}
