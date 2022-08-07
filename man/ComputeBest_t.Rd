\name{ComputeBest_t}
\alias{ComputeBest_t}
\title{
  Monte Carlo simulation to investigate the optimal number of
  points to use in the moment conditions
}
\description{
  Runs Monte Carlo simulation for different values of \eqn{\alpha} and
  \eqn{\beta} and computes a specified number of t-points that minimises
  the determinant of the asymptotic covariance matrix.
}
\usage{
ComputeBest_t(AlphaBetaMatrix = abMat, nb_ts = seq(10, 100, 10),
              alphaReg = 0.001, FastOptim = TRUE, ...)
}
\arguments{
  \item{AlphaBetaMatrix}{
    values of the parameter \eqn{\alpha} and \eqn{\beta} from which we
    simulate the data. By default, the values of \eqn{\gamma} and
    \eqn{\delta} are set to 1 and 0, respectively;
    a \eqn{2 \times  n} matrix.
  }
  \item{nb_ts}{
    vector of numbers of t-points to use for the minimisation;
    default = \code{seq(10, 100, 10)}.
  }
  \item{alphaReg}{
    value of the regularisation parameter; numeric, default = 0.001.
  }
  \item{FastOptim}{
    Logical flag; if set to TRUE, \code{optim} with "Nelder-Mead" method
    is used (fast but not accurate). Otherwise, \code{nlminb} is used
    (more accurate but slower).
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
  \code{\link{ComputeBest_tau}},
  \code{\link{Best_t-class}}
}
\keyword{Simulation}
