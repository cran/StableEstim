\name{RegularisedSol}
\alias{RegularisedSol}
\title{
  Regularised Inverse
}
\description{

  Regularised solution of the (ill-posed) problem \eqn{K\phi = r} where
  \eqn{K} is a \eqn{n \times n} matrix, \eqn{r} is a given vector of
  \code{length} n. Users can choose one of the 3 schemes described in
  Carrasco and Florens (2007).

}
\usage{
RegularisedSol(Kn, alphaReg, r,
               regularization = c("Tikhonov", "LF", "cut-off"),
               ...)
}
\arguments{
  \item{Kn}{
    numeric \eqn{n \times n} matrix.
  }
  \item{alphaReg}{
    regularisation parameter; numeric in ]0,1].
  }
  \item{r}{
    numeric vector of \code{length} n.
  }
  \item{regularization}{
    regularization scheme to be used, one of \code{"Tikhonov"} (Tikhonov
    scheme), \code{"LF"} (Landweber-Fridmann) and \code{"cut-off"}
    (spectral cut-off). See Details.
  }
  \item{\dots}{
    the value of \eqn{c} used in the \code{"LF"} scheme. See
    Carrasco and Florens(2007).
  }
}
\details{

  Following Carrasco and Florens(2007), the regularised solution of the
  problem \eqn{K \phi=r} is given by :
  \deqn{\varphi_{\alpha_{reg}} =
    \sum_{j=1}^{n} q(\alpha_{reg},\mu_j)\frac{<r,\psi_j >}{\mu_j} \phi_j
    ,
  }
  where \eqn{q} is a (positive) real function with some regularity
  conditions and \eqn{\mu,\phi,\psi} the singular decomposition of the
  matrix \eqn{K}.

  The \code{regularization} parameter defines the form of the function
  \eqn{q}. For example, the \code{"Tikhonov"} scheme defines
  \eqn{q(\alpha_{reg},\mu) = \frac{\mu^2}{\alpha_{reg}+\mu^2}}.

  When the matrix \eqn{K} is symmetric, the singular decomposition is
  replaced by a spectral decomposition.
}
\value{
  the regularised solution, a vector of length n.
}
\references{ % bibentry: carrasco2007linear
Carrasco M, Florens J and Renault E (2007).
``Linear inverse problems in structural econometrics estimation based on spectral decomposition and regularization.''
\emph{Handbook of econometrics}, \bold{6}, pp. 5633--5751.
% end:bibentry:  carrasco2007linear
}
\seealso{
\code{\link{solve}}
}
\examples{
## Adapted from R examples for Solve 
## We compare the result of the regularized sol to the expected solution

hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+")}

K_h8 <- hilbert(8);
r8 <- 1:8

alphaReg_robust <- 1e-4
Sa8_robust <- RegularisedSol(K_h8,alphaReg_robust,r8,"LF")

alphaReg_accurate <- 1e-10
Sa8_accurate <- RegularisedSol(K_h8,alphaReg_accurate,r8,"LF")

## when pre multiplied by K_h8, the expected solution is 1:8
## User can check the influence of the choice of alphaReg
}
\keyword{general-functions}
