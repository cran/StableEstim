\name{TexSummary}
\alias{TexSummary}
\title{
  LaTex summary.
}
\description{
  Creates a Tex table from a summary object or a vector of files.
}
\usage{
TexSummary(obj, files = NULL, sep_ = ",", FctsToApply = StatFcts,
           caption = "Statistical Summary",label = "Simtab",
           digits = 3, par_index = 1,MCparam=1000,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
    \code{list} of \code{length} 4 containing a summary \code{matrix} object
    associated to each parameter identical to the one produced by the
    function \code{\link{ComputeStatObjectFromFiles}}.
  }
  \item{files}{
    \code{character} vector containing the files name to be parsed. It
    will be passed to function \code{\link{ComputeStatObjectFromFiles}}.
  }
  \item{sep_}{
    Field separator character to be passed to the function
    \code{\link{ComputeStatObjectFromFiles}}.
  }
  \item{FctsToApply}{
    Functions used to produce the statistical summary to be passed
    to the function \code{\link{ComputeStatObjectFromFiles}}.
  }
  \item{caption}{
    \code{character} vector with length equal to length \code{par_index}
    containing the table's caption or title.
  }
  \item{label}{
    \code{Character} vector with length equal to length \code{par_index}
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
    \code{numeric} or \code{character} vector of length 1,2,3 or 4 of the
    desired indexes to be selected in \code{obj}. See details}
  \item{MCparam}{ Number of Monte Carlo simulation for each couple of parameter,
      default=1000; integer
    }
  \item{\dots}{
    Other arguments to be passed to the function
    \code{\link{ComputeStatObjectFromFiles}}.
  }
}
\details{
  Accepted values for \code{par_index} are \code{c(1,2,3,4)} or
  \code{c("alpha","beta","gamma","delta")} or mixed.
  %\cr\cr
  Some examples are provided in the example folder.
  %\cr\cr
}
\value{
  Returns a \code{list} of objects (same length as \code{par_index}) from class
  \code{Latex}. See  \code{toLatex}  
}
\seealso{
  \code{\link{Estim_Simulation}},
  \code{\link{ComputeStatObjectFromFiles}},
  \code{xtable}
}

\keyword{Simulation}
