\name{sampleComplexCFMoment}
\alias{sampleComplexCFMoment}

\title{
  Complex moment condition based on the characteristic function.
}
\description{
  Computes the moment condition based on the characteristic function as a complexl vector.
}
\usage{
sampleComplexCFMoment(x, t, theta, pm = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{
   vector of data where the ecf is computed.
 }
 \item{t}{
   Vector of (real) numbers where the CF is evaluated; numeric
 }
 \item{theta}{	
   Vector of parameters of the stable law; vector of length 4.
 }	  
  \item{pm}{
    Parametrisation, an integer (0 or 1); default: \code{pm=0}(
    the Nolan \sQuote{S0} parametrisation).
  }
}
\details{
  \bold{The moment conditions}
  %\cr\cr
  The moment conditions are given by:
  \deqn{g_t(X,\theta)=g(t,X;\theta)= e^{itX} - \phi_{\theta}(t)}
  If one has a sample \eqn{x_1,\dots,x_n} of i.i.d realisations of the same random variable \eqn{X}, then:
  \deqn{\hat{g}_n(t,\theta)  = \frac{1}{n}\sum_{i=1}^n g(t,x_i;\theta) =  \phi_n(t) -\phi_\theta(t)}
  where \eqn{\phi_n(t)} is the eCF associated to the sample
  \eqn{x_1,\dots,x_n} and defined by \eqn{\phi_n(t)= \frac{1}{n} \sum_{j=1}^n e^{itX_j}}
  %\cr\cr
  The function compute the vector of difference between the eCF and the CF at a set of given point \code{t}.
}
\value{
  Returns a complex vector of \code{length(t)}.
}
\seealso{
  	\code{\link{ComplexCF}},\code{\link{sampleRealCFMoment}}
}
\examples{
## define the parameters
 nt <- 10   
 t <- seq(0.1,3,length.out=nt)
 theta <- c(1.5,0.5,1,0)
pm <- 0

set.seed(222);x=rstable(200,theta[1],theta[2],theta[3],theta[4],pm)

# Compute the characteristic function
CFMC <- sampleComplexCFMoment(x=x,t=t,theta=theta,pm=pm)

}

\keyword{stable-functions}
