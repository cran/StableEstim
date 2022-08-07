\name{TexSummary}
\alias{TexSummary}
\title{
  LaTeX summary
}
\description{
  Creates a TeX table from a summary object or a vector of files.
}
\usage{
TexSummary(obj, files = NULL, sep_ = ",", FctsToApply = StatFcts,
           caption = "Statistical Summary", label = "Simtab",
           digits = 3, par_index = 1, MCparam = 1000, ...)
}
\arguments{
  \item{obj}{
    \code{list} of length 4 containing a summary \code{matrix}
    object associated to each parameter identical to the one produced by
    function \code{\link{ComputeStatObjectFromFiles}}.
  }
  \item{files}{
    \code{character} vector containing the files name to be parsed. It
    will be passed to function \code{\link{ComputeStatObjectFromFiles}}.
  }
  \item{sep_}{
    field separator character to be passed to function
    \code{\link{ComputeStatObjectFromFiles}}.
  }
  \item{FctsToApply}{
    functions used to produce the statistical summary to be passed
    to the function \code{\link{ComputeStatObjectFromFiles}}.
  }
  \item{caption}{
    \code{character} vector with length equal to \code{length(par_index)}
    containing the table's caption or title.
  }
  \item{label}{
    \code{character} vector with length equal to \code{length(par_index)} 
    containing the LaTeX label.
  }
  \item{digits}{
    \code{numeric} vector of length equal to one (in which case it will
    be replicated as necessary) or to the number of columns of the
    resulting table or length of \code{FctsToApply} or matrix of the
    same size as the resulting table indicating the number of digits
    to display in the corresponding columns. See \code{xtable}.
  } 
  \item{par_index}{
    \code{numeric} or \code{character} vector of length 1, 2, 3 or 4 of
    the desired indices to be selected in \code{obj}. See Details.
  }
  
  \item{MCparam}{
    number of Monte Carlo simulations for each couple of parameters,
    default = 1000; integer.
  }
  \item{\dots}{
    other arguments to be passed to function
    \code{\link{ComputeStatObjectFromFiles}}.
  }
}
\details{

  Accepted values for \code{par_index} are \code{c(1,2,3,4)} or
  \code{c("alpha","beta","gamma","delta")} or mixed.

  Some examples are provided in the example folder.

}
\value{
  a \code{list} of length \code{length(par_index)} whose elements are
  objects from class \code{Latex} (produced by \code{toLatex})
}
\seealso{
  \code{\link{Estim_Simulation}},
  \code{\link{ComputeStatObjectFromFiles}},
  \code{xtable}
}
\keyword{Simulation}
