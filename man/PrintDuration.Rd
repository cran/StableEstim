\name{PrintDuration}
\alias{PrintDuration}
\title{
 Print duration
}
\description{
  Print duration in human readable format.
}
\usage{
PrintDuration(t, CallingFct = "")
}
\arguments{
  \item{t}{
    Duration; \code{numeric} of \code{length} 1 or 3.
  }
  \item{CallingFct}{
    Name of the calling function.
  }
}
\details{
  The duration will be printed in the format: hours/minutes/seconds.
}
\value{
  Prints a \code{character} to the screen.
}

\examples{
ti <- getTime_()
for (i in 1:100) x = i*22.1
tf <- getTime_()
duration <- ComputeDuration(ti, tf)
PrintDuration(duration, "test")
}
\keyword{general-functions}
