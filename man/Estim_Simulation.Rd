\name{Estim_Simulation}
\alias{Estim_Simulation}
\title{
Monte Carlo Simulation
}
\description{
Runs Monte Carlo simulation for a selected estimation method. The
function can save a file and produce a statistical summary.
}
\usage{
Estim_Simulation(AlphaBetaMatrix = abMat, SampleSizes = c(200, 1600),
                 MCparam = 100, Estimfct = c("ML", "GMM", "Cgmm","Kout"),
                 HandleError = TRUE, FctsToApply = StatFcts,
                 saveOutput = TRUE, StatSummary = FALSE,
                 CheckMat=TRUE,tolFailCheck=tolFailure,
                 SeedOptions=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AlphaBetaMatrix}{
    Values of the parameter \eqn{\alpha} and  \eqn{\beta} from which
    we simulate the data. By default, The values of \eqn{\gamma}
    and \eqn{\delta} are set to 1 and 0 respectively; matrix \eqn{2 \times n}
  }
  \item{SampleSizes}{
    Sample sizes to be used to simulate the data. By default, we use
    \code{200} (small sample size) and \code{1600} (large sample size);
      vector of integer.
    }
    \item{MCparam}{
      Number of Monte Carlo simulation for each couple of parameter,
      default=100; integer
    }
    \item{Estimfct}{
      The estimation function to be used. User can choose between
      \code{ML}, \code{GMM}, \code{Cgmm} or \code{Kout}.
    }
    \item{HandleError}{
      Logical flag: if set to TRUE, the simulation doesn't stop when an
      error in the estimation function is encountered. A vector of (size 4)
      \code{NA} is saved and the the simulation carries on. See details.
    }
    \item{FctsToApply}{
      Functions used to produce the statistical summary. See details ;
      vector of \code{character}.
    }
    \item{saveOutput}{
      Logical flag: if set to TRUE, a csv file (for each couple of
      parameter \eqn{\alpha} and \eqn{\beta}) with the the estimation
      information is saved in the current directory. See details
    }
    \item{StatSummary}{
      Logical flag: if set to TRUE, a statistical summary (using
      \code{FctsToApply}) is returned. See details.
    }
    \item{CheckMat}{
      Logical flag: if set to TRUE, an estimation is
      declared failed if the squared error of the estimation is larger
      than \code{tolFailCheck}; default TRUE
    }
    \item{tolFailCheck}{
      Tolerance on the squared error of the estimation
      to be declared failed; default 1.5
    }
    \item{SeedOptions}{
      List to control the seed generation. See details.
      }
    \item{\dots}{
      Other arguments to be passed to the estimation function.
    }
  }
\details{
  %\cr\cr
  \bold{Error Handling}
  %\cr\cr
  It is advisable to set it to TRUE when user is planning to launch long
  simulations as it will prevent the procedure to stop if an error
  occurs for one sample data. The estimation function will produce a vector of
  \code{NA} as estimated parameters related to this (error
  generating) sample data and move on to the next Monte Carlo step.
  %\cr\cr
  \bold{Statistical summary}
  %\cr\cr
  The function is able to produce a statistical summary of the Monte
  Carlo simulation for each parameter (slices of the list). Each slice
  is a matrix where the rows represents the true values of the
  parameters and the columns the statistical information.
  %\cr\cr
  In all cases, the following quantities are computed:
  \describe{
    \item{\code{sample size}:}{ the sample size used to produce the simulated data}
    \item{\code{alphaT}, \code{betaT}:}{ the true values of the
      parameters}
    \item{\code{failure}:}{ the number of times the procedure failed to
      produce relevant estimation.}
    \item{\code{time}:}{ The average running time in seconds of the
      estimation procedure}
  }
  Besides, the (vector of \code{character}) \code{FctsToApply} controls
  the other quantities to be computed by providing the name of the
  function object to be applied to the vector of estimated
  parameter. The signature of the function should be of the form
  \code{fctName=function(p,...){...}} where \code{p} is the vector
  (\code{length(p)=MCparam}) of parameter estimate and \code{\dots} is
  the extra arguments to be passed the function.
  %\cr\cr
  By default, the functions of the \code{StatFcts} will be
  applied but user can pass his own functions by providing their name (
  in the \code{FctsToApply} vector and their definition in the global environment.
  %\cr\cr
  Note that if \code{CheckMat} is set to TRUE, the estimation is
  considered failed if the squared error (of the first 2 parameters
  \code{alpha} and \code{beta}) is larger than \code{tolFailCheck}. 
  \bold{Output file}
  %\cr\cr
  Setting \code{saveOutput} to TRUE will have the side effect of
  saving a csv file in the working directory. This file will have
  \code{MCparam}*\code{length(SampleSizes)} lines and its columns will
  be :
  \describe{
    \item{\code{alphaT}, \code{betaT}:}{ the true value of the parameters.}
    \item{\code{data size}:}{ the sample size used to generate the simulated data.}
    \item{\code{seed}:}{ the seed value used to generate the simulated data.}
    \item{\code{alphaE}, \code{betaE}, \code{gammaE}, \code{deltaE}:}{ the
    estimate of the 4 parameters.}
    \item{\code{failure}:}{binary: 0 for success, 1 for failure.}
    \item{\code{time}:}{ estimation running time in seconds.}
  }
  The file name is informative to let the user identify the value
  of the true parameters, the MC parameters as well as the options
  selected for the estimation method.
  %\cr\cr
  The csv file is updated after each MC estimation which is useful when
  the simulation stops before it finishes. Besides, using the
  check-pointing mechanism explained below, the simulation can re-start
  from where it stopped.
  %\cr\cr
  \emph{Check-pointing}
  Checkpointing is the act of saving enough program state and results
  so far calculated that a computation can be stopped and restarted. The
  way we did it here is to save a text file with some useful
  information about the state of the estimation. This text file is
  updated after each MC iteration and read at the beginning of function
  \code{Estim_Simulation} to allow the simulation to re-start from where
  it stopped. This file is deleted at the end of the simulation procedure.
  %\Cr\cr
  \emph{SeedOptions}
  If user does not want to control the seed generation, he could ignore
  this argument (default value \code{NULL}). This argument can be more
  useful when one wants to cut the simulation (even for one parameter
  value) into pieces. In that case, he can control which part of the
  seed vector he wants to use.
  \describe{
    \item{\code{MCtot}:}{ total values of MC simulations in the entire process.}
    \item{\code{seedStart}:}{ starting index in the seed vector. The
      vector extracted will be of size \code{MCparam}.}
  }
}
\value{
  If \code{StatSummary} is set to TRUE, a \code{list} with 4 slices
  (corresponding to the 4 parameters) is returned. Each slice is a
  matrix. If \code{SaveOutput} is set to TRUE, only a csv file
  is saved and nothing is returned (if \code{StatSummary} is FALSE).
  If both are FALSE, the function stops.
}

\seealso{
  \code{\link{Estim}},\code{\link{CgmmParametersEstim}},
  \code{\link{GMMParametersEstim}},\code{\link{MLParametersEstim}}
}


\keyword{Simulation}

