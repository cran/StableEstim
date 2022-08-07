\name{Estim-class}
\Rdversion{1.1}
\docType{class}
\alias{Estim-class}
\alias{initialize,Estim-method}
\alias{show,Estim-method}

\title{Class \code{"Estim"}}
\description{

  Class for storing the results of estimating parameters of stable
  laws, output of function \code{Estim()}.

}
\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("Estim", par,
  ...)}.  Users can provide some (or all) of the inputs stated below to
  create an object from this class or call function \code{\link{Estim}}
  with appropriate arguments.

}
\section{Slots}{
  \describe{
    \item{\code{par}:}{
      \code{numeric(4)}, values of the 4 estimated parameters.
    }
    \item{\code{par0}:}{
      \code{numeric(4)}, initial values for the 4 parameters.      
    }
    \item{\code{vcov}:}{
      object of class \code{"matrix"} (\code{4 x 4}),
      representing the covariance matrix of the estimated parameters.
    }
    \item{\code{confint}:}{
      object of class \code{"matrix"} (\code{4 x 4}),
      representing the confidence interval computed
      at a specific level (attribute of the object).
    }
    \item{\code{data}:}{
      \code{numeric()}, the data used to compute the estimates.
    }
    \item{\code{sampleSize}:}{
      \code{numeric(1)}, length of the data.
    }
    \item{\code{others}:}{
      \code{list()},
      further information about the estimation method.
    }
    \item{\code{duration}:}{
      \code{numeric(1)}, duration in seconds.
    }
    \item{\code{failure}:}{
      \code{numeric(1)},
      represents the status of the procedure: 0 failure or 1 success.
    }
    \item{\code{method}:}{
      Object of class \code{"character"},
      description of the parameters used in the estimation.
    }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "Estim")}:
      creates an object of this class using the inputs described above.
    }
    \item{show}{\code{signature(object = "Estim")}:
      summarised print of the object.
    }
  }
}
\seealso{
  \code{\link{Estim}}
}
\keyword{classes}
