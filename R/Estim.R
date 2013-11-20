Estim <- function(EstimMethod=c("ML","GMM","Cgmm","Kout"),data,theta0=NULL,
                  ComputeCov=FALSE,HandleError=TRUE,...){
    
    if (is.null(theta0)) theta0 <- IGParametersEstim(data,...)
    OutputObj <- new(Class="Estim",par0=theta0,data=data,failure=1)
    if (missing(data)) stop("data not provided !")
    method <- match.arg(EstimMethod)

    EstimFcts <- getEstimFcts(method)
    

    res <- .initRes(method)    
    if(HandleError){
        tr <- tryCatch(EstimFcts$Params(x=data,theta0=theta0,...),
                       error=function(e)e) 
        err <- inherits(tr, "error")
        if (!err) {res <- tr ;OutputObj@failure <- 0}
    }
    else {res <- EstimFcts$Params(x=data,theta0=theta0,...);OutputObj@failure <- 0}
    OutputObj@par <- NameParamsObjects(res$Estim$par)
    OutputObj@others <- res$Estim
    OutputObj@duration <- as.numeric(res$duration)
    OutputObj@method <- res$method
    

    if (ComputeCov) {
        OutputObj@vcov <- EstimFcts$CovarianceMat(data=OutputObj@data,
                                                  EstimObj=res,
                                                  ...)
                                                  
        OutputObj@confint <- AsymptoticConfidenceInterval(thetaEst=OutputObj@par,
                                                          n_sample=OutputObj@sampleSize,
                                                          Cov=OutputObj@vcov,
                                                          qLaw=qnorm,...) 
    }
    
    OutputObj
}

Estim_Des<- function(EstimMethod=c("ML","GMM","Cgmm","Kout"),...){
    method <- match.arg(EstimMethod)

    EstimFcts <- getEstimFcts(method)    
    EstimFcts$methodDes(...)
}

getEstimFcts <- function(method=c("ML","GMM","Cgmm","Kout")){
    Output <- switch(method,
                     "ML"={list(Params=MLParametersEstim,
                                CovarianceMat=.asymptoticVarianceEstimML,
                                methodDes=.methodDesML
                                )
                       },
                     "GMM"={list(Params=GMMParametersEstim,
                                 CovarianceMat=.asymptoticVarianceEstimGMM,
                                 methodDes= getGMMmethodName
                                 )
                        },
                     "Cgmm"={list(Params=CgmmParametersEstim,
                                  CovarianceMat=.asymptoticVarianceEstimCgmm,
                                  methodDes=getCgmmMethodName
                                  )
                         },
                     "Kout"={list(Params=KoutParametersEstim,
                                 CovarianceMat=.asymptoticVarianceEstimKout,
                                 methodDes=getKoutMethodName
                                 )
                     },
                     stop(paste(method," not taken into account !"))
                     )
    Output
}

.initRes <- function(method){
    list(Estim=list(par= rep(NaN,4)),
         duration=0,
         method=paste(method,"failed",sep="_")
         )
}

## Matrix: 4x2 (col1=min,col2=max)
AsymptoticConfidenceInterval <- function(thetaEst,n_sample,Cov,qLaw=qnorm,level=0.95,...){
    mat <- matrix(NaN,ncol=2,nrow=4);attr(mat,"level") <- level
    z   <- qLaw(level) 
    
    V   <- sqrt(diag(Cov))*z
    mat[,1] <- thetaEst-V
    mat[,2] <- thetaEst+V
    
    NameParamsObjects(mat)
}
