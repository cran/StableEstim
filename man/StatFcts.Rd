\name{StatFcts}
\alias{StatFcts}
\docType{data}
\title{
  Default functions used to produce the statistical summary
}
\description{
  Default functions used to produce the statistical summary in the
  Monte Carlo simulations.
}
\details{
  The functions are:
  \describe{
    \item{Mean}{\code{.mean <- function(p,...) mean(p)}}
    \item{Min}{\code{.min <- function(p,...) min(p)}}
    \item{Max}{\code{.max <- function(p,...) max(p)}}
    \item{Sn}{\code{.Sn <- function(p,n,...) sqrt(n)*sd(p)}}
    \item{MSE}{\code{.MSE <- function(p,paramT,...) (1/length(p))*sum((p-paramT)^2) }}
    \item{Std error}{\code{.st.err <- function(p,...) sd(p)/sqrt(length(p))}}
  }
  To change the statistical summary, provide functions with similar
  signatures and pass a character vector containing the function
  names to \code{\link{Estim_Simulation}}.
}
\keyword{general-functions}
