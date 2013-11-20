\name{ConcatFiles}
\alias{ConcatFiles}
\title{
  Concatenates output files.
}
\description{
  Creates a unique file by concatenating several output files associated
  to one set of parameters.
}
\usage{
ConcatFiles(files, sep_ = ",", outfile, headers_ = TRUE,
            DeleteIfExists=TRUE)
}

\arguments{
  \item{files}{
    \code{character} Vector containing the files name to be concatenated.
    See details.
  }
  \item{sep_}{
    Field separator character to be used in function
    \code{read.csv()} and \code{write.table()}.
    Values on each line of the file are separated by
    this character; It can also be a vector character
    (same length as \code{files}) if different separators are useed for
    each file; default: ","
  }
  \item{outfile}{
    Name of the output file; \code{character}
  }
  \item{headers_}{
    Vector of \code{boolean} of length 1 or same length as \code{files}
    to indicate for each file if the header argument is to be considered
    or not. To be passed to function \code{read.csv()}. 
  }
  \item{DeleteIfExists}{if \code{outfile} exists, it will be deleted and
    recreated (over-written).
  }
}
\details{
  The files to be concatenated should be related to the same set of
  parameters \code{alpha} and \code{beta}. The function stops if one of
  the file contains 2 (or more) different set of parameters (the
  function compares the values of columns 1 and 2 row by row) or if the
  set of parameters within one file is different from the one from other
  files.
}

\value{
  Returns an output file \code{outfile} saved in the working directory.
}
\seealso{
   \code{\link{Estim_Simulation}}
}

\keyword{Simulation}
