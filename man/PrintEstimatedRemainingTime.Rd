\name{PrintEstimatedRemainingTime}
\alias{PrintEstimatedRemainingTime}
\title{
  Estimated remaining time
}
\description{
  Prints the estimated remaining time in a loop.
  Useful in Monte Carlo simulations.
}
\usage{
PrintEstimatedRemainingTime(ActualIter, ActualIterStartTime, TotalIterNbr)
}
\arguments{
  \item{ActualIter}{
    Actual Iteration; \code{integer}
  }
  \item{ActualIterStartTime}{
    Actual Iteration Starting time; \code{numeric}
  }
  \item{TotalIterNbr}{
    Total number of iterations; \code{integer}
  }	
}
\details{
  Called at the end of each Monte Carlo step, this function will compute
  the duration of the actual step, an estimate of the remaining MC loops
  duration and prints the result to the screen in a human readable
  format using function \code{\link{PrintDuration}}.
}
\seealso{
  \code{\link{PrintDuration}},
  \code{\link{ComputeDuration}}.
}
\keyword{general-functions}
