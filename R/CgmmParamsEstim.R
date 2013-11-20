CgmmParametersEstim <- function(x,type=c("2S","IT","Cue"),alphaReg=0.01,
                                subdivisions=50,IntegrationMethod=c("Uniform","Simpson"),
                                randomIntegrationLaw=c("unif","norm"),s_min=0,s_max=1,
                                theta0=NULL,IterationControl=list(),pm=0,
                                PrintTime=FALSE,...){
    
    if (is.null(theta0)) theta0 <- IGParametersEstim(x,pm) 
    type <- match.arg(type)
    
    t_init <- getTime_()
    method <- getCgmmMethodName(type=type,alphaReg=alphaReg,
                                subdivisions=subdivisions,
                                IntegrationMethod=IntegrationMethod,
                                randomIntegrationLaw=randomIntegrationLaw,
                                s_min=s_min,s_max=s_max
                            )
    
    Estim <- switch(type,
                    "2S"={Compute2SCgmmParametersEstim(x=x,theta0=theta0,alphaReg=alphaReg,
                                                       pm=pm,s_min=s_min,s_max=s_max,
                                                       IntegrationMethod=IntegrationMethod,
                                                       randomIntegrationLaw=randomIntegrationLaw,
                                                       subdivisions=subdivisions,...)},
                    "IT"={ComputeITCgmmParametersEstim(x=x,theta0=theta0,alphaReg=alphaReg,
                                                       pm=pm,s_min=s_min,s_max=s_max,
                                                       IntegrationMethod=IntegrationMethod,
                                                       randomIntegrationLaw=randomIntegrationLaw,
                                                       subdivisions=subdivisions,
                                                       IterationControl=IterationControl,...)},
                    "Cue"={ComputeCueCgmmParametersEstim(x=x,theta0=theta0,alphaReg=alphaReg,
                                                         pm=pm,s_min=s_min,s_max=s_max,
                                                         IntegrationMethod=IntegrationMethod,
                                                         randomIntegrationLaw=randomIntegrationLaw,
                                                         subdivisions=subdivisions,
                                                         IterationControl=IterationControl,...)},                    
                    stop(paste(type," not taken into account for Cgmm procedure"))
                    )
    if (PrintTime){
        CallingFct <- paste("CgmmParametersEstim",type,sep="_")
        PrintDuration(ComputeDuration(t_init,getTime_()),CallingFct)
    }
    list(Estim=Estim,
         duration=as.numeric(ComputeDuration(t_init,getTime_(),TRUE)),
         method=method)
}

Compute2SCgmmParametersEstim <- function(x,theta0,alphaReg,pm,s_min,
                                         s_max,IntegrationMethod,randomIntegrationLaw,
                                         subdivisions,...){

    dots <- list(...)
    if (is.null(dots$control)){
        control <- list(abs.tol=1e-15,
                        rel.tol=1e-7,
                        x.tol=1.5e-5,
                        xf.tol=2.2e-10)
        thetaHat <- as.numeric(nlminb(start=theta0,objective=ComputeObjectiveCgmm,
                                      gradient=NULL,hessian=NULL,
                                      Weighting="Id",x=x,
                                      alphaReg=alphaReg,pm=pm,
                                      thetaHat=NULL,s_min=s_min,s_max=s_max,
                                      IntegrationMethod=IntegrationMethod,
                                      randomIntegrationLaw=randomIntegrationLaw,
                                      subdivisions=subdivisions,...,
                                      control=control,
                                      lower = c(eps, -1+eps, 0+eps, -Inf),
                                      upper = c(2-eps, 1-eps,  Inf,  Inf))$par)
    }
    else{
        thetaHat <- as.numeric(nlminb(start=theta0,objective=ComputeObjectiveCgmm,
                                      gradient=NULL,hessian=NULL,
                                      Weighting="Id",x=x,
                                      alphaReg=alphaReg,pm=pm,
                                      thetaHat=NULL,s_min=s_min,s_max=s_max,
                                      IntegrationMethod=IntegrationMethod,
                                      randomIntegrationLaw=randomIntegrationLaw,
                                      subdivisions=subdivisions,...,
                                      lower = c(eps, -1+eps, 0+eps, -Inf),
                                      upper = c(2-eps, 1-eps,  Inf,  Inf))$par)
    }
    
    Cmat <- ComputeCmat(x=x,thetaHat=thetaHat,s_min=s_min,s_max=s_max,
                        IntegrationMethod=IntegrationMethod,
                        randomIntegrationLaw=randomIntegrationLaw,
                        subdivisions=subdivisions,pm=pm)

    if (is.null(dots$control)){
        control <- list(abs.tol=1e-15,
                        rel.tol=1e-7,
                        x.tol=1.5e-5,
                        xf.tol=2.2e-10)
        
        res <- nlminb(start=theta0,objective=ComputeObjectiveCgmm,
                      Weighting="optimal",Cmat=Cmat,x=x,
                      alphaReg=alphaReg,pm=pm,
                      thetaHat=thetaHat,s_min=s_min,s_max=s_max,
                      IntegrationMethod=IntegrationMethod,
                      randomIntegrationLaw=randomIntegrationLaw,
                      subdivisions=subdivisions,...,
                      control=control,
                      lower = c(eps, -1+eps, 0+eps, -Inf),
                      upper = c(2-eps, 1-eps,  Inf,  Inf))
    }
    else{
        res <- nlminb(start=theta0,objective=ComputeObjectiveCgmm,
                      Weighting="optimal",Cmat=Cmat,x=x,
                      alphaReg=alphaReg,pm=pm,
                      thetaHat=thetaHat,s_min=s_min,s_max=s_max,
                      IntegrationMethod=IntegrationMethod,
                      randomIntegrationLaw=randomIntegrationLaw,
                      subdivisions=subdivisions,...,
                      lower = c(eps, -1+eps, 0+eps, -Inf),
                      upper = c(2-eps, 1-eps,  Inf,  Inf))
    }
    list(par=as.numeric(res$par),all=res)
}

ComputeITCgmmParametersEstim <- function(x,theta0,alphaReg,pm,s_min,
                                         s_max,IntegrationMethod,randomIntegrationLaw,
                                         subdivisions,IterationControl,...){

    iter=0
    IterationControl <- checkIterationControl(IterationControl)
    
    theta1 <- as.numeric(nlminb(start=theta0, objective=ComputeObjectiveCgmm,
                                Weighting="Id",x=x,
                                alphaReg=alphaReg,pm=pm,
                                thetaHat=NULL,s_min=s_min,s_max=s_max,
                                IntegrationMethod=IntegrationMethod,
                                randomIntegrationLaw=randomIntegrationLaw,
                                subdivisions=subdivisions,...,
                                lower = c(eps, -1+eps, 0+eps, -Inf),
                                upper = c(2-eps, 1-eps,  Inf,  Inf))$par)
    
    PrevEstimParVal <- theta1
    RelativeErr=IterationControl$RelativeErrMax+5
    
    while((iter < IterationControl$NbIter) && (RelativeErr > IterationControl$RelativeErrMax) ){
        Cmat <- ComputeCmat(x=x,thetaHat=PrevEstimParVal,s_min=s_min,s_max=s_max,
                            IntegrationMethod=IntegrationMethod,
                            randomIntegrationLaw=randomIntegrationLaw,
                            subdivisions=subdivisions,pm=pm)
        dots <- list(...)
        if (is.null(dots$control)){
            control <- list(abs.tol=1e-15,
                            rel.tol=1e-7,
                            x.tol=1.5e-5,
                            xf.tol=2.2e-10)
            CurrentEstimAllInfo <- nlminb(start=PrevEstimParVal,
                                          objective=ComputeObjectiveCgmm, 
                                          Weighting="optimal",Cmat=Cmat,
                                          alphaReg=alphaReg,pm=pm,x=x,
                                          thetaHat=PrevEstimParVal,s_min=s_min,s_max=s_max,
                                          IntegrationMethod=IntegrationMethod,
                                          randomIntegrationLaw=randomIntegrationLaw,
                                          subdivisions=subdivisions,...,
                                          control=control,
                                          lower = c(eps, -1+eps, 0+eps, -Inf),
                                          upper = c(2-eps, 1-eps,  Inf,  Inf))
        }
        else{
            CurrentEstimAllInfo <- nlminb(start=PrevEstimParVal,
                                          objective=ComputeObjectiveCgmm, 
                                          Weighting="optimal",Cmat=Cmat,
                                          alphaReg=alphaReg,pm=pm,x=x,
                                          thetaHat=PrevEstimParVal,s_min=s_min,s_max=s_max,
                                          IntegrationMethod=IntegrationMethod,
                                          randomIntegrationLaw=randomIntegrationLaw,
                                          subdivisions=subdivisions,...,
                                          lower = c(eps, -1+eps, 0+eps, -Inf),
                                          upper = c(2-eps, 1-eps,  Inf,  Inf))
        }
        CurrentEstimParVal <- CurrentEstimAllInfo$par 
        
        if (IterationControl$PrintIter) PrintIteration(CurrentEstimParVal,iter,IterationControl$NbIter) 

        RelativeErr <- abs(CurrentEstimParVal-PrevEstimParVal)
        PrevEstimParVal <- CurrentEstimParVal
        iter=iter+1
    }
    
    list(par=as.numeric(CurrentEstimParVal),all=CurrentEstimAllInfo)
}

ComputeCueCgmmParametersEstim <- function(x,theta0,alphaReg,pm,s_min,
                                         s_max,IntegrationMethod,randomIntegrationLaw,
                                         subdivisions,IterationControl,...){

    iter=0
    IterationControl <- checkIterationControl(IterationControl)
    
    theta1 <- as.numeric(nlminb(start=theta0, objective=ComputeObjectiveCgmm,
                                Weighting="Id",x=x,
                                alphaReg=alphaReg,pm=pm,
                                thetaHat=NULL,s_min=s_min,s_max=s_max,
                                IntegrationMethod=IntegrationMethod,
                                randomIntegrationLaw=randomIntegrationLaw,
                                subdivisions=subdivisions,...,
                                lower = c(eps, -1+eps, 0+eps, -Inf),
                                upper = c(2-eps, 1-eps,  Inf,  Inf))$par)
    
    PrevEstimParVal <- theta1
    RelativeErr=IterationControl$RelativeErrMax+5
    
    while((iter < IterationControl$NbIter) && (RelativeErr > IterationControl$RelativeErrMax) ){
        dots <- list(...)
        if (is.null(dots$control)){
            control <- list(abs.tol=1e-15,
                            rel.tol=1e-7,
                            x.tol=1.5e-5,
                            xf.tol=2.2e-10)
            
            CurrentEstimAllInfo <- nlminb(start=PrevEstimParVal,
                                          objective=ComputeObjectiveCgmm,
                                          Cmat=NULL,Weighting="optimal",
                                          alphaReg=alphaReg,pm=pm,x=x,
                                          thetaHat=PrevEstimParVal,s_min=s_min,s_max=s_max,
                                          IntegrationMethod=IntegrationMethod,
                                          randomIntegrationLaw=randomIntegrationLaw,
                                          subdivisions=subdivisions,...,
                                          control=control,
                                          lower = c(eps, -1+eps, 0+eps, -Inf),
                                          upper = c(2-eps, 1-eps,  Inf,  Inf))
        }
        else{
            CurrentEstimAllInfo <- nlminb(start=PrevEstimParVal,
                                          objective=ComputeObjectiveCgmm,
                                          Cmat=NULL,Weighting="optimal",
                                          alphaReg=alphaReg,pm=pm,x=x,
                                          thetaHat=PrevEstimParVal,s_min=s_min,s_max=s_max,
                                          IntegrationMethod=IntegrationMethod,
                                          randomIntegrationLaw=randomIntegrationLaw,
                                          subdivisions=subdivisions,...,
                                          lower = c(eps, -1+eps, 0+eps, -Inf),
                                          upper = c(2-eps, 1-eps,  Inf,  Inf))
        }
        CurrentEstimParVal <- as.numeric(CurrentEstimAllInfo$par)
        
        if (IterationControl$PrintIter) PrintIteration(CurrentEstimParVal,iter,IterationControl$NbIter) 

        RelativeErr <- abs(CurrentEstimParVal-PrevEstimParVal)
        PrevEstimParVal <- CurrentEstimParVal
        iter=iter+1
    }
    
    list(par=as.numeric(CurrentEstimParVal),all=CurrentEstimAllInfo)
}

ComputeObjectiveCgmm <- function(theta,Cmat=NULL,x,Weighting=c("optimal","Id"),
                                 alphaReg,pm,thetaHat,s_min,s_max,
                                 subdivisions=50,IntegrationMethod=c("Uniform","Simpson"),
                                 randomIntegrationLaw=c("norm","unif"),...){
    
    n <- length(x)
    IntegrationMethod <- match.arg(IntegrationMethod)
    randomIntegrationLaw <- match.arg(randomIntegrationLaw)
    Weighting <- match.arg(Weighting)

    ObjectiveVal <-ComputeCgmmFcts(Fct="Objective",theta=theta,Cmat=Cmat,x=x,
                                   Weighting=Weighting,alphaReg=alphaReg,pm=pm,
                                   thetaHat=thetaHat,s_min=s_min,s_max=s_max,
                                   subdivisions=subdivisions,
                                   IntegrationMethod=IntegrationMethod,
                                   randomIntegrationLaw=randomIntegrationLaw,...)
    as.numeric(Mod(ObjectiveVal))
}

## we mean the inverse of the CovarianceCgmm called V in the article 
ComputeCovarianceCgmm <- function(theta,Cmat=NULL,x,alphaReg,pm,thetaHat,s_min,s_max,
                                  subdivisions=50,IntegrationMethod=c("Uniform","Simpson"),
                                  randomIntegrationLaw=c("norm","unif"),...){
    
    n <- length(x)
    IntegrationMethod <- match.arg(IntegrationMethod)
    randomIntegrationLaw <- match.arg(randomIntegrationLaw)

    CovMat <- ComputeCgmmFcts(Fct="Covariance",theta=theta,Cmat=Cmat,x=x,
                              Weighting="optimal",alphaReg=alphaReg,pm=pm,
                              thetaHat=thetaHat,s_min=s_min,s_max=s_max,
                              subdivisions=subdivisions,
                              IntegrationMethod=IntegrationMethod,
                              randomIntegrationLaw=randomIntegrationLaw,...)
    CovMat/(n-4)
}

.asymptoticVarianceEstimCgmm <- function(data,EstimObj,...){
    V <- ComputeCovarianceCgmm(theta=EstimObj$Estim$par,
                               thetaHat=EstimObj$Estim$par,
                               x=data,...)
    NameParamsObjects(Mod(ComputeCutOffInverse(V))/length(data))
}



ComputeCgmmFcts <- function(Fct=c("Objective","Covariance"),theta,Cmat=NULL,x,
                            Weighting=c("optimal","Id"),
                            alphaReg,pm,thetaHat,s_min,s_max,
                            subdivisions=50,IntegrationMethod=c("Uniform","Simpson"),
                            randomIntegrationLaw=c("norm","unif"),...){
    
    n <- length(x)
    Fct <- match.arg(Fct)
    IntegrationMethod <- match.arg(IntegrationMethod)
    randomIntegrationLaw <- match.arg(randomIntegrationLaw)
    Weighting <- match.arg(Weighting)
    
    ghatBarFctOft <-function(t_var,X) Conj(sampleComplexCFMoment(x=X,t=t_var,theta=theta,pm=pm))
    ghatFctOft <-function(t_var,X) Conj(ghatBarFctOft(t_var,X))
    
    if(Weighting=="Id"){
        ObjectiveVal <- IntegrateRandomVectorsProduct(f_fct=ghatFctOft,X=x,g_fct=ghatBarFctOft,
                                                      Y=x,s_min=s_min,s_max=s_max,
                                                      subdivisions=subdivisions,
                                                      IntegrationMethod=IntegrationMethod,
                                                      randomIntegrationLaw=randomIntegrationLaw,...)
    }
    else{
        V <- ComputeV(Fct=Fct,theta=theta,thetaHat=thetaHat,X=x,
                      IntegrationMethod=IntegrationMethod,s_min=s_min,
                      randomIntegrationLaw=randomIntegrationLaw,s_max=s_max,
                      subdivisions=subdivisions,pm=pm,...)
        if (is.null(Cmat)){
            if(Fct=="Covariance") thetaToUse <- thetaHat
            else thetaToUse <- theta
            Cmat <- ComputeCmat(x=x,thetaHat=thetaToUse,s_min=s_min,s_max=s_max,IntegrationMethod=IntegrationMethod,
                                randomIntegrationLaw=randomIntegrationLaw,subdivisions=subdivisions,pm=pm,...)/(n-4)
        }
        else Cmat <- Cmat/(n-4)

        In <- diag(nrow=n,ncol=n)
        Cmat2 <- Cmat %*%Cmat
        matrixToInverse <- alphaReg*In + Cmat2
        ObjectiveVal <- crossprod(Conj(V),solve(a=matrixToInverse,b=V))
    }
    ObjectiveVal
}

ComputeCmat <- function(x,thetaHat,s_min,s_max,IntegrationMethod,randomIntegrationLaw,subdivisions,pm,...){
    f_fct <- function(s,x){sapply(X=x,FUN=sampleComplexCFMoment,t=s,theta=thetaHat,pm=pm)}
    f_bar_fct <- function(s,x){Conj(f_fct(s,x))}
    
    IntegrateRandomVectorsProduct(f_fct=f_bar_fct,X=x,g_fct=f_fct,Y=x,s_min=s_min,s_max=s_max,
                                  subdivisions=subdivisions,IntegrationMethod=IntegrationMethod,
                                  randomIntegrationLaw=randomIntegrationLaw,...)
}

ComputeV <- function(Fct=c("Objective","Covariance"),theta,thetaHat,X,
                     s_min,s_max,IntegrationMethod,
                     randomIntegrationLaw,subdivisions,
                     pm,...){

    Fct <- match.arg(Fct)
    
    g_hat_fct <- function(s,x){sampleComplexCFMoment(x=x,t=s,theta=theta,pm=pm)}
    g_bar_fct <- function(s,x){Conj(sapply(X=x,FUN=sampleComplexCFMoment,t=s,theta=thetaHat,pm=pm))}
    Jac_g_hat_fct <- function(s,x){jacobianSampleComplexCFMoment(t=s,theta=theta,pm=pm)}

    if (Fct=="Covariance"){
        res <- IntegrateRandomVectorsProduct(f_fct=g_bar_fct,X=X,g_fct=Jac_g_hat_fct,Y=X,s_min=s_min,s_max=s_max,
                                             subdivisions=subdivisions,IntegrationMethod=IntegrationMethod,
                                             randomIntegrationLaw=randomIntegrationLaw,...)
    }
    else if (Fct=="Objective"){
        res <- IntegrateRandomVectorsProduct(f_fct=g_bar_fct,X=X,g_fct=g_hat_fct,Y=X,s_min=s_min,s_max=s_max,
                                             subdivisions=subdivisions,IntegrationMethod=IntegrationMethod,
                                             randomIntegrationLaw=randomIntegrationLaw,...)
    }
    res
}
 
getCgmmMethodName <- function(type,alphaReg,
                              subdivisions,IntegrationMethod,
                              randomIntegrationLaw,s_min,s_max,...){
    args <- list(...)
    paste("Cgmm",
          paste("type=",type,sep=""),
          paste("alphaReg=",alphaReg,sep=""),
          paste("OptimAlgo=","nlminb",sep=""),
          paste("subdivisions=",subdivisions,sep=""),
          paste("IntegrationMethod=",IntegrationMethod,sep=""),
          paste("randomIntegrationLaw=",randomIntegrationLaw,sep=""),
          paste("s_min=",s_min,sep=""),
          paste("s_max=",s_max,sep=""),
          sep="_")    
}


