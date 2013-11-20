MLParametersEstim<-function(x,theta0=NULL,pm=0,PrintTime=FALSE,...)
{
    t_init<- getTime_()
    if (is.null(theta0)) theta0 <- IGParametersEstim(x,pm)
    method <- .methodDesML()

    dots <- list(...)
    if (is.null(dots$control)){
        control <- list(factr=1e5,
                        pgtol=1e-5)
        
        Estim <- optim(par=theta0,fn=SumLogDesnity,
                       gr=NULL,x=x,pm=pm,
                       control=control,...,
                       method="L-BFGS-B",
                       lower = c(eps, -1+eps, 0+eps, -Inf),
                       upper = c(2-eps, 1-eps,  Inf,  Inf)
                       )
    }
    else{
        Estim <- optim(par=theta0,fn=SumLogDesnity,
                       gr=NULL,x=x,pm=pm,...,
                       method="L-BFGS-B",
                       lower = c(eps, -1+eps, 0+eps, -Inf),
                       upper = c(2-eps, 1-eps,  Inf,  Inf)
                       )
    }
    
    if (PrintTime) PrintDuration(ComputeDuration(t_init,
                                                 t_final <- getTime_()
                                                 ),
                                 "MLParametersEstim")
    
    list(Estim=Estim,duration=ComputeDuration(t_init,getTime_(),TRUE),method=method)
}

.methodDesML <- function(...){
    l <- list(...)
    paste("ML",
          paste("OptimAlgo=","L-BFGS-B",sep=""),
          sep="_")
}

SumLogDesnity <- function(theta,x,sign=-1,pm=0,...)
    {
        sign*sum(dstable(x,alpha = theta[1], beta = theta[2],
                         gamma = theta[3], delta = theta[4],
                         pm=pm,log=TRUE))
    }

VectorialDensity <- function(theta,xi){
    dstable(xi,theta[1],theta[2],
            theta[3],theta[4])
}

## the result is a mxn matrix i.e lenght(x) x 4
jacVectorialDensity <- function(theta,xi){
    NumDeriv_jacobian(fctToDeriv=VectorialDensity,
                      WhereFctIsEvaluated=theta,
                      xi=xi)
}

.plot.integrand <- function(theta,xmin=-20,xmax=20,type=1){
    x <- seq(from=xmin,to=xmax,by=0.2)
    invf <- 1/VectorialDensity(theta,x)
    df <- jacVectorialDensity(theta,x)
    y <- invf*df[,type]*df[,type]
    plot(x,y,type="l")
}

## asymptotic covariance C=1/n (I)^-1
## I_ij=int_{-\infty}^{\infty}df/d\theta_i df/d\theta_j 1/f dx
## the function is SLOW as we need to compute the
## integral accuractly using the QUADPACK routines
## dqags and dqagi called by fct integrate

invFisherMatrix <- function(theta,subdivisions=100){
    mat <- matrix(NA,4,4)
    integrand <- function(x,i,j)
        {
            invf <- 1/VectorialDensity(theta,x)
            df <- jacVectorialDensity(theta,x)
            y <- invf*df[,i]*df[,j]
        }
    
    for (i in 1:4){
        for (j in 1:i){
            mat[i,j] <- integrate(f=integrand,lower=-Inf,upper=Inf,
                   i=i,j=j,subdivisions=subdivisions)$value
            mat[j,i] <- mat[i,j]
        }
    }
    solve(mat)
}

asymptoticVarianceEstimML <- function(thetaEst,n_sample,subdivisions=100,...){
       NameParamsObjects(invFisherMatrix(as.numeric(thetaEst),subdivisions)/n_sample)
   }

.asymptoticVarianceEstimML <- function(data,EstimObj,...){
   asymptoticVarianceEstimML(thetaEst=EstimObj$Estim$par,
                             n_sample=length(data),
                             ...)
}
