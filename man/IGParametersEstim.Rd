\name{IGParametersEstim}
\alias{IGParametersEstim}
\title{
  Estimate parameters of stable laws by Kogon and McCulloch methods
}
\description{
  Kogon regression method is used together with the McCulloch quantile
  method to provide initial estimates of parameters of stable
  distributions. 
}
\usage{
IGParametersEstim(x, pm = 0, ...)
}

\arguments{
  \item{x}{
    data used to perform the estimation: vector of length n.
  }
  \item{pm}{
    parametrisation, an integer (0 or 1); default: \code{pm = 0}
    (Nolan's \sQuote{S0} parametrisation).
  }
  \item{\dots}{
    other arguments. Currently not used.
  }
}
\details{
  The parameters \eqn{\gamma} and \eqn{\delta} are estimated using the
  McCulloch(1986) quantile method from \pkg{fBasics}. The data is
  rescaled using those estimates and used to perform the Kogon
  regression method to estimate \eqn{\alpha} and \eqn{\beta}.
}
\value{
  a vector of length 4 containing the estimates of the 4 parameters.
}
\references{
  % bibentry: kogon1998characteristic
Kogon SM and Williams DB (1998).
``Characteristic function based estimation of stable distribution parameters.''
\emph{A practical guide to heavy tailed data}, pp. 311--335.
% end:bibentry:  kogon1998characteristic
  % bibentry: mcculloch1986simple
McCulloch JH (1986).
``Simple consistent estimators of stable distribution parameters.''
\emph{Communications in Statistics-Simulation and Computation}, \bold{15}(4), pp. 1109--1136.
% end:bibentry:  mcculloch1986simple
}
\seealso{
  \code{\link{Estim}}, \code{\link{McCullochParametersEstim}}
}
\examples{
x <- rstable(200, 1.2, 0.5, 1, 0, pm = 0)
IGParametersEstim(x, pm = 0)
}
\keyword{Estim-functions}
