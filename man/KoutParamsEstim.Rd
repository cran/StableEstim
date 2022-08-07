\name{KoutParametersEstim}
\alias{KoutParametersEstim}
\title{
  Iterative Koutrouvelis regression method
}
\description{
  Iterative Koutrouvelis regression method with different
  spacing schemes (points where the eCF is computed).
}
\usage{
KoutParametersEstim(x, theta0 = NULL,
                    spacing = c("Kout", "UniformSpac", "ArithSpac", "free"),
                    pm = 0, tol = 0.05, NbIter = 10, PrintTime = FALSE, ...)
}
\arguments{
  \item{x}{
    data used to perform the estimation: vector of length n.
  }
  \item{theta0}{
    initial guess for the 4 parameters values: vector of length 4
  }
  \item{spacing}{
    scheme used to select the points where the moment conditions are
    evaluated. \code{Kout} is the scheme suggested by Koutrouvelis,
    \code{UniformSpac} and \code{ArithSpac} are the uniform and
    arithmetic spacing schemes over the informative interval
    [\eqn{\epsilon},\eqn{A_n}]. If user choose free, he needs to provide
    a set of points \code{t_points} and \code{u_points} in \code{\dots}.
  }
  \item{pm}{
    parametrisation, an integer (0 or 1); default: \code{pm = 0}
    (Nolan's \sQuote{S0} parametrisation).
  }
  \item{tol}{
    the loop stops if the relative error between two consecutive
    estimation is smaller then \code{tol}; default = 0.05.
  }
  \item{NbIter}{
    maximum number of iteration. The loop stops when \code{NbIter}
    is reached; default = 10.
  }
  \item{PrintTime}{
    logical flag; if set to TRUE, the estimation duration is printed out
    to the screen in a readable format (h/min/sec).
  }
  \item{\dots}{
    other arguments to pass to the function. See Details.
  }
}
\details{
  \bold{spacing}

  4 options for the spacing scheme are implemented as described
  above. In particular:
  \describe{
    \item{\code{UniformSpac}, \code{ArithSpac}:}{
      The user can specify the number of points to choose in both
      regression by inputting \code{nb_t} and \code{nb_u}. Otherwise the
      Koutrouvelis table will be used to compte them.
    }
    \item{\code{free}:}{
      The user is expected to provide \code{t_points} and
      \code{u_points} otherwise the \code{Kout} scheme will be used.
    }
  }
}
\value{
  a list with the following elements:
  \item{Estim}{
    \code{list} containing the vector of 4 parameters estimate
    (\code{par}), the 2 regressions objects (\code{reg1} and
    \code{reg2}) and the matrix of iterations estimate (\code{vals}).
  }
  \item{duration}{estimation duration in a numerical format.}
  \item{method}{\code{character} describing the method used.}
}
\references{% bibentry: koutrouvelis1980regression
Koutrouvelis IA (1980).
``Regression-type estimation of the parameters of stable laws.''
\emph{Journal of the American Statistical Association}, \bold{75}(372), pp. 918--928.
% end:bibentry:  koutrouvelis1980regression

   % bibentry: koutrouvelis1981iterative
Koutrouvelis IA (1981).
``An iterative procedure for the estimation of the parameters of stable laws: An iterative procedure for the estimation.''
\emph{Communications in Statistics-Simulation and Computation}, \bold{10}(1), pp. 17--28.
% end:bibentry:  koutrouvelis1981iterative
}
\seealso{
  \code{\link{Estim}}
}
\examples{
pm <- 0
theta <- c(1.45, 0.5, 1.1, 0.4)
set.seed(1235)
x <- rstable(200, theta[1], theta[2], theta[3], theta[4], pm = pm)
theta0 <- theta - 0.1
spacing <- "Kout"

KoutParametersEstim(x = x, theta0 = theta0, 
                        spacing = spacing, pm = pm)
}
\keyword{Estim-functions}
