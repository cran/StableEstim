\name{Estim-class}
\Rdversion{1.1}
\docType{class}
\alias{Estim-class}
\alias{initialize,Estim-method}
\alias{show,Estim-method}

\title{Class \code{"Estim"}}
\description{
  Class storing all the information about the estimation method; output
  of function Estim 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Estim", par, ...)}.
User can provide some/all of the inputs stated below to create an object
of this class or call function \code{\link{Estim}} with appropriate arguments. 
}
\section{Slots}{
  \describe{
    \item{\code{par}:}{Object of class \code{"numeric"}; Value of the
      4 estimated parameters.}
    \item{\code{par0}:}{Object of class \code{"numeric"}; Initial guess
      for the 4 parameters. }
    \item{\code{vcov}:}{Object of class \code{"matrix"} with 4 rows and
      4 columns representing the covariance matrix }
    \item{\code{confint}:}{Object of class \code{"matrix"} with 4 rows
      and 2 columns representing the confidence interval computed
      at a specific level (attribute of the object).}
    \item{\code{data}:}{Object of class \code{"numeric"} used to compute
    the estimation.}
    \item{\code{sampleSize}:}{Object of class \code{"numeric"} ; length
      of the data }
    \item{\code{others}:}{Object of class \code{"list"} ; more
      information about the estimation method }
    \item{\code{duration}:}{Object of class \code{"numeric"} ; duration
      in seconds }
    \item{\code{failure}:}{Object of class \code{"numeric"} representing
    the status of the procedure: 0 failure or 1 success}
    \item{\code{method}:}{Object of class \code{"character"} description
    of the parameter used in the estimation.}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "Estim")}: creates un
      aboject of this class using the inputs described above }
    \item{show}{\code{signature(object = "Estim")}: summarised print of
      the object. }
	 }
}

\seealso{
  \code{\link{Estim}}
}

\keyword{classes}
