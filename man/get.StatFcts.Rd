\name{get.StatFcts}
\alias{get.StatFcts}
\title{
  Default functions used to produce the statistical summary
}
\description{
  Default functions used to produce the statistical summary in the
  Monte Carlo simulations.
}
\usage{
get.StatFcts()
}
\value{
  The functions computed are:
  \describe{
    \item{Mean}{\code{.mean <- function(p,...) mean(p)}}
    \item{Min}{\code{.min <- function(p,...) min(p)}}
    \item{Max}{\code{.max <- function(p,...) max(p)}}
    \item{Sn}{\code{.Sn <- function(p,n,...) sqrt(n)*sd(p)}}
    \item{MSE}{\code{.MSE <- function(p,paramT,...) (1/length(p))*sum((p-paramT)^2)}}
    \item{Std error}{\code{.st.err <- function(p,...) sd(p)/sqrt(length(p))}}
  }
  Users can define their own summaries by defining functions with
  similar signatures and passing a \code{character} vector containing
  the functions' names to \code{\link{Estim_Simulation}}.
}
\keyword{data-functions}
