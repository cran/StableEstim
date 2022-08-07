\name{ComputeDuration}
\alias{ComputeDuration}
\title{
  Duration
}
\description{
  Compute the duration between 2 time points.
}
\usage{
ComputeDuration(t_init, t_final, OneNumber = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t_init}{Starting time; numeric.}
  \item{t_final}{Final time; numeric.}
  \item{OneNumber}{
    Logical flag; if set to TRUE, the duration in seconds will be
    returned. Otherwise, a vector of length 3 will be computed
    representing the time in h/min/sec.
  }
}
\value{
  a \code{numeric} of length 1 or 3 depending on the value of
  \code{OneNumber} flag.
}
\seealso{
  \code{\link{PrintDuration}},
  \code{\link{PrintEstimatedRemainingTime}}.
}
\examples{
ti <- getTime_()
for (i in 1:100) x <- i*22.1
tf <- getTime_()
ComputeDuration(ti,tf)
}
\keyword{general-functions}
