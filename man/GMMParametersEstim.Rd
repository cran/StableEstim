\name{GMMParametersEstim}
\alias{GMMParametersEstim}
\title{
generalised method of moment with finite number of moment conditions
}
\description{
Generalised method of moment with finite number of moment conditions. It
uses a regularisation technique to make the method more robust (when
the number of moment condition is large) and allows different schemes to
select where the moment conditions are computed.
}
\usage{
GMMParametersEstim(x, algo = c("2SGMM", "ITGMM", "CueGMM"),
                   alphaReg =0.01,
                   regularization = c("Tikhonov", "LF", "cut-off"),
                   WeightingMatrix = c("OptAsym", "DataVar", "Id"),
                   t_scheme = c("equally", "NonOptAr", "uniformOpt",
                                "ArithOpt", "VarOpt", "free"),
                   theta0 = NULL,
                   IterationControl = list(),
                   pm = 0, PrintTime = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Data used to perform the estimation: vector of length n.
  }
  \item{algo}{
    GMM algorithm: "2SGMM" is the two step GMM proposed by Hansen(1982)
    and the "CueGMM" and "ITGMM" are respectively the continuous updated and the
    iterative GMM proposed by Hansen, Eaton et Yaron (1996) and adapted
    to the continuum case.
  }
  \item{alphaReg}{
    Value of the regularisation parameter; numeric,default=0.01
  }
  \item{regularization}{
    Regularization scheme to be used. User can choose between
    \code{"Tikhonov"} scheme ,\code{"LF"} (Landweber-Fridmann) and
    \code{"cut-off"}(spectral cut-off). See \code{\link{RegularisedSol}}.
  }
  \item{WeightingMatrix}{
    Type of Weighting matrix used to compute the objective function. User
    can choose between \code{"OptAsym"} (the optimal asymptotic),
    \code{"DataVar"} (the data driven) and \code{"Id"} (the identity
    matrix). See details
  }
  \item{t_scheme}{
    Scheme used to select the points where the moment conditions are
    evaluated. User can choose between \code{"equally"} (equally placed),
    \code{"NonOptAr"} (non optimal arithmetic placement),
    \code{"uniformOpt"} (uniform optimal placement), \code{"ArithOpt"}
    (arithmetic optimal placement) , \code{"Var Opt"}(optimal variance
    placement) and \code{"free"} (user needs to pass his own set of point
    in the \code{\dots}). See details.
  }
  \item{theta0}{
    Initial guess for the 4 parameters values: If \code{NULL}, the
    Kogon-McCulloch method is called, see
    \code{\link{IGParametersEstim}}; vector of length 4.
  }
  \item{IterationControl}{
    Only used if \code{type="IT"} or \code{type="Cue"} to control the
    iterations. See details
  }
  \item{pm}{
    Parametrisation, an integer (0 or 1); default: \code{pm=0} (the
    Nolan \sQuote{S0} parametrisation).
  }
  \item{PrintTime}{
    Logical flag; if set to TRUE, the estimation duration is printed out to
    the screen in a readable format (h/min/sec).
  }
  \item{\dots}{Other arguments to pass to the regularisation function,
    the optimisation function or the selection scheme (including the
    function that finds the first zero of the eCF). See details.
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
  \bold{Objective function}
  %\cr\cr
  \deqn{obj{\theta}= < K^{-1/2} \hat{g}_n(.;\theta),K^{-1/2}\hat{g}_n(.;\theta)>}
  where \eqn{K^{-1}f} denotes the solution \eqn{\varphi} (when it
  exists) of the equation \eqn{K \varphi=f} and \eqn{K^{-1/2}=(K^{-1})^{1/2}}.
  The optimal choice of the Weighting operator K (a matrix in the GMM case) and
  its estimation are discussed in Hansen (1982).
  %\cr\cr
  \bold{Weighting operator (Matrix)}
  %\cr\cr
  \describe{
    \item{\code{OptAsym}:}{the optimal asymptotic choice as described by
    Hansen. The expression of the components of this matrix could be
    found for example in Feuerverger and McDunnough (1981b)}
    \item{\code{DataVar}:}{the covariance matrix of the data provided}
    \item{\code{Id}:}{the identity matrix}
  }
  %\cr\cr
  \bold{the t-scheme}
  %\cr\cr
  One of the most improtant feature of this method is that it allows the
  user to choose how to place the points where the moment conditions are
  evaluated. The general rule is that user can provide his own set of
  points (option \code{"free"}) or choose one of the other scheme. In the
  latter case he \emph{needs to specify the number of points}
  \code{nb_t} in the \code{\dots} argument and eventually the lower and
  upper limit (by setting \code{Constrained} to FALSE and
  providing \code{min_t} and \code{max_t}) in the non-optimised case. If
  one of the optimised case is selected, setting \code{Constrained} to FALSE
  will not constrain the choice of \eqn{\tau}, see below.
  We mean by optimised set of point, the set that minimises the
  (determinant) of the asymptotic covariance matrix as suggested by
  Schmidt (1982) and Besbeas and Morgan (2008).
  %\cr
  6 options have been implemented:
  \describe{
    \item{\code{"equally"}:}{equally placed point in
      [\code{min_t},\code{max_t}]. When provided, user's \code{min_t}
      and \code{max_t} will be used (when\code{Coinstrained}=FALSE).
      Otherwise, \code{eps} and \code{An} will be used instead (where
      \code{An} is the first zero of the eCF)}
    \item{\code{"NonOptAr"}:}{non optimal arithmetic
      placement: \eqn{t_j=\frac{j(j+1)}{nbt(nbt+1)} (max-eps);j=1,\dots,nbt} where max
    is the upper band of the set of points selected as discussed before.}
    \item{\code{"uniformOpt"}:}{uniform optimal placement: \eqn{t_j=j \tau, j=1,\dots nbt}}
    \item{\code{"ArithOpt"}:}{arithmetic optimal placement: \eqn{t_j=j(j+1) \tau, j=1,\dots nbt}}
    \item{\code{"Var Opt"}:}{optimal variance placement as explained above.}
    \item{\code{"free"}:}{user needs to pass his own set of point in the \code{\dots}}
  }
  For the \code{"ArithOpt"} and \code{"uniformOpt"} schemes, the function
  to minimise is seen as a function of the real parameter \eqn{\tau}
  instead of doing a vectorial optimisition as in the \code{"Var Opt"}
  case. In the latter case, one can choose between a fast (but less
  accurate) optimisation routine or a slow (but more accurate) one by
  setting the \code{FastOptim} flag to the desired value.
  %\cr\cr
  \bold{The IterationControl}
  %\cr\cr
  If \code{type="IT"} or \code{type="Cue"} user can control each
  iteration by setting up the \code{list} \code{IterationControl} which
  contains the following elements:
  \describe{
    \item{\code{NbIter}:}{maximum number of iteration. The loop stops when \code{NBIter}
      is reached; default=10}
    \item{\code{PrintIterlogical}:}{if set to TRUE, the value of the current
      parameter estimation is printed to the screen at each iteration; default=TRUE}
    \item{\code{RelativeErrMax}:}{the loop stops if the relative error between two
      consecutive estimation is smaller than
      \code{RelativeErrMax}; default=1e-3}
  }
}

\value{
  Returns a list with the following elements:
  \item{Estim}{output of the optimisation function}
  \item{duration}{estimation duration in a numerical format}
  \item{method}{\code{character} describing the method used}
  \item{tEstim}{final set of points selected for the
    estimation. Only relevant when one of the optimisation scheme is selected.}
}
\references{
  % bibentry: hansen1982large
Hansen LP (1982).
``Large sample properties of generalized method of moments estimators.''
\emph{Econometrica: Journal of the Econometric Society}, pp. 1029--1054.
% end:bibentry:  hansen1982large

  % bibentry: hansen1996finite
Hansen LP, Heaton J and Yaron A (1996).
``Finite-sample properties of some alternative GMM estimators.''
\emph{Journal of Business \& Economic Statistics}, \bold{14}(3), pp. 262--280.
% end:bibentry:  hansen1996finite

  % bibentry: feuerverger1981efficient
Feuerverger A and McDunnough P (1981).
``On efficient inference in symmetric stable laws and processes.''
\emph{Statistics and Related Topics}, \bold{99}, pp. 109--112.
% end:bibentry:  feuerverger1981efficient

  % bibentry: feuerverger1981some
Feuerverger A and McDunnough P (1981).
``On some Fourier methods for inference.''
\emph{Journal of the American Statistical Association}, \bold{76}(374), pp. 379--387.
% end:bibentry:  feuerverger1981some

  % bibentry: schmidt1982improved
Schmidt P (1982).
``An improved version of the Quandt-Ramsey MGF estimator for mixtures of normal distributions and switching regressions.''
\emph{Econometrica: Journal of the Econometric Society}, pp. 501--516.
% end:bibentry:  schmidt1982improved

  % bibentry: besbeas2008improved
Besbeas P and Morgan B (2008).
``Improved estimation of the stable laws.''
\emph{Statistics and Computing}, \bold{18}(2), pp. 219--231.
% end:bibentry:  besbeas2008improved

}

\note{\code{nlminb} was used for the minimisation of the GMM objective
  funcion and to compute \eqn{tau} in the \code{"uniformOpt"} and
  \code{"ArithOpt"} schemes. In the \code{"Var Opt"} scheme,
  \code{optim} was preferred. All those routines have been
  selected after running different tests using the summary table produced
  by package \pkg{optimx} for comparing the performance of different
  optimisation methods.
}

\seealso{\code{\link{Estim}}}

\examples{
## General data
theta <- c(1.5,0.5,1,0)
pm <- 0
set.seed(345);
x <- rstable(100,theta[1],theta[2],theta[3],theta[4],pm)
##---------------- 2S free ----------------
## method specific arguments
regularization="cut-off"
WeightingMatrix="OptAsym"
alphaReg=0.005
t_seq=seq(0.1,2,length.out=12)

## If you are just interested by the value
## of the 4 estimated parameters
t_scheme="free"
algo = "2SGMM"
     
suppressWarnings(GMMParametersEstim(x=x,
                                    algo=algo,alphaReg=alphaReg,
                                    regularization=regularization,
                                    WeightingMatrix=WeightingMatrix,
                                    t_scheme=t_scheme,
                                    pm=pm,PrintTime=TRUE,t_free=t_seq))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Estim-functions}

