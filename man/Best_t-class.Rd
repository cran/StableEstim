\name{Best_t-class}
\Rdversion{1.1}
\docType{class}
\alias{Best_t-class}
\alias{+,Best_t,Best_t-method}
\alias{initialize,Best_t-method}
\alias{show,Best_t-method}

\title{Class \code{"Best_t"}}
\description{
  Class used to store the result of function \code{\link{ComputeBest_t}}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("Best_t", theta, nbt, tvec, detVal, convcode, ...)},
  where the user can specify some/all of the inputs or call function
  \code{\link{ComputeBest_t}}.
}
\section{Slots}{
  \describe{
    \item{\code{theta}:}{Object of class \code{"vector"}; values of the
      4 parameters.}
    \item{\code{nbt}:}{Object of class \code{"vector"}; number of points
      used in the minimisation.}
    \item{\code{tvec}:}{Object of class \code{"list"}; values of the
      best t-vectors.}
    \item{\code{detVal}:}{Object of class \code{"vector"}; values of the
      optimal determinant found after minimisation.}
    \item{\code{convcode}:}{Convergence code.}
  }
}
\section{Methods}{
  \describe{
    \item{+}{\code{signature(e1 = "Best_t", e2 = "Best_t")}:
      sum objects from class \code{Best_t}.
    }
    \item{initialize}{\code{signature(.Object = "Best_t")}:
      initialise an object from class \code{Best_t} as described above.
    }
    \item{show}{\code{signature(object = "Best_t")}:
      print a summary of the object.
    }
  }
}
\seealso{
  \code{\link{ComputeBest_t}}
}
\keyword{classes}
