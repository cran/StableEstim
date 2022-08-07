\name{McCullochParametersEstim}
\alias{McCullochParametersEstim}
\title{
  Quantile-based method
}
\description{
  McCulloch quantile-based method.
}
\usage{
McCullochParametersEstim(x)
}
\arguments{
  \item{x}{
    data used to perform the estimation: vector of length n.
  }
}
\details{
  The code is a modified version of function \code{.qStableFit} from
  package \pkg{fBasics}. 
}
\value{
  \code{numeric} of length 4, represening the value of the 4
  parameters. 
}
\references{
  % bibentry: mcculloch1986simple
McCulloch JH (1986).
``Simple consistent estimators of stable distribution parameters.''
\emph{Communications in Statistics-Simulation and Computation}, \bold{15}(4), pp. 1109--1136.
% end:bibentry:  mcculloch1986simple
}
\seealso{
  \code{\link{Estim}},
  \code{\link{IGParametersEstim}}
}
\examples{
set.seed(333)
x <- rstable(500, 1.3, 0.4, 1, 0)
McCullochParametersEstim(x)
}
\keyword{Estim-functions}
