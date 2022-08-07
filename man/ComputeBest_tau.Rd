\name{ComputeBest_tau}
\alias{ComputeBest_tau}
\title{
  Run Monte Carlo simulation to investigate the optimal \eqn{\tau}
}
\description{
  Runs Monte Carlo simulation to investigate the optimal number of
  points to use when one of the reduced spacing schemes is considered.
}
\usage{
ComputeBest_tau(AlphaBetaMatrix = abMat, nb_ts = seq(10, 100, 10),
                tScheme = c("uniformOpt", "ArithOpt"),
                Constrained = TRUE, alphaReg = 0.001, ...)
}
\arguments{
  \item{AlphaBetaMatrix}{
    values of the parameter \eqn{\alpha} and \eqn{\beta} from which we
    simulate the data. By default, the values of \eqn{\gamma} and
    \eqn{\delta} are set to 1 and 0, respectively;
    a \eqn{2 \times n} matrix.
  }
  \item{nb_ts}{
    vector of number of t-points to use for the minimisation;
    default = \code{seq(10,100,10)}.
  }
  \item{tScheme}{
    scheme used to select the points where the moment conditions are
    evaluated, one of \code{"uniformOpt"} (uniform optimal placement)
    and \code{"ArithOpt"} (arithmetic optimal placement). See function
    \code{\link{GMMParametersEstim}}.
  }
  \item{Constrained}{
    logical flag: if set to True, lower and upper bands will be computed
    as discussed for function \code{\link{GMMParametersEstim}}.
  }
  \item{alphaReg}{
    value of the regularisation parameter; numeric, default = 0.001.
  }
  \item{\dots}{
    Other arguments to pass to the optimisation function.
  }
}
\value{
  a \code{list} containing slots from class \code{\link{Best_t-class}}
  corresponding to one value of the parameters \eqn{\alpha} and
  \eqn{\beta}.
}
\seealso{
  \code{\link{ComputeBest_t}},
  \code{\link{Best_t-class}}
}
\keyword{Simulation}
