\name{ComputeBest_tau}
\alias{ComputeBest_tau}

\title{
  Runs Monte Carlo Simulation to investigate the optimal \eqn{\tau}.
}
\description{
  Runs Monte Carlo Simulation to investigate the optimal number of
  points to use when one of the reduced spacing scheme is considered.
}
\usage{
ComputeBest_tau(AlphaBetaMatrix = abMat, nb_ts = seq(10, 100, 10),
                tScheme = c("uniformOpt", "ArithOpt"),
                Constrained = TRUE, alphaReg = 0.001, ...)
}

\arguments{
  \item{AlphaBetaMatrix}{
    Values of the parameter \eqn{\alpha} and  \eqn{\beta} from which
    we simulate the data. By default, The values of \eqn{\gamma}
    and \eqn{\delta} are set to 1 and 0 respectively; matrix \eqn{2 \times n}
  }
  \item{nb_ts}{
    Vector of number of t-points to use for the minimisation;
    default \code{seq(10,100,10)}.
  }
  \item{tScheme}{
    Scheme used to select the points where the moment conditions are
    evaluated. User can choose between \code{"uniformOpt"}
    (uniform optimal placement) and  \code{"ArithOpt"} (arithmetic optimal
    placement). See function \code{\link{GMMParametersEstim}}
  }
  \item{Constrained}{
    Logical flag: if set to True lower and upper bands will be computed as
    discussed function \code{\link{GMMParametersEstim}}.
  }
  \item{alphaReg}{
    Value of the regularisation parameter; numeric, default=0.001.
  }
  \item{\dots}{
    Other arguments to pass to the optimisation function.
  }
}

\value{
  Returns a \code{list} containing slots from class \code{\link{Best_t-class}}
  corresponding to one value of the parameters \eqn{\alpha} and
  \eqn{\beta}.
}
\seealso{
  \code{\link{ComputeBest_t}},\code{\link{Best_t-class}}
}

\keyword{Simulation}% __ONLY ONE__ keyword per line
