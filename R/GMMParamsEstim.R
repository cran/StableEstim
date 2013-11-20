GMMParametersEstim <- function(x,algo=c("2SGMM","ITGMM","CueGMM"),
                               alphaReg=0.01,regularization=c("Tikhonov","LF","cut-off"),
                               WeightingMatrix=c("OptAsym","DataVar","Id"),
                               t_scheme=c("equally","NonOptAr","uniformOpt","ArithOpt","VarOpt","free"),
                               theta0=NULL,IterationControl=list(),pm=0,
                               PrintTime=FALSE,...){

    if (is.null(theta0)) theta0 <- IGParametersEstim(x,pm) 
    algo <- match.arg(algo)
    regularization <- match.arg(regularization)
    WeightingMatrix <- match.arg(WeightingMatrix)
    t_scheme <- match.arg(t_scheme)

    t_init <- getTime_()
    method <- getGMMmethodName(algo=algo,alphaReg=alphaReg,
                               regularization=regularization,
                               WeightingMatrix=WeightingMatrix,
                               t_scheme=t_scheme,...)
    
    Estim <- switch(algo,
                    "2SGMM"={Compute2SGMMParametersEstim(x=x,theta0=theta0,alphaReg=alphaReg,
                                                         regularization=regularization,
                                                         WeightingMatrix=WeightingMatrix,
                                                         t_scheme=t_scheme,pm=pm,...)},
                    "ITGMM"={ComputeITGMMParametersEstim(x=x,theta0=theta0,alphaReg=alphaReg,
                                                         regularization=regularization,
                                                         WeightingMatrix=WeightingMatrix,
                                                         t_scheme=t_scheme,IterationControl=IterationControl,
                                                         pm=pm,...)},
                    "CueGMM"={ComputeCueGMMParametersEstim(x=x,theta0=theta0,alphaReg=alphaReg,
                                                           regularization=regularization,
                                                           WeightingMatrix=WeightingMatrix,
                                                           t_scheme=t_scheme,IterationControl=IterationControl,
                                                           pm=pm,...)},
                    stop(paste(algo," not taken into account !"))
                    )
    if (PrintTime){
        CallingFct <- paste("GMMParametersEstim",algo,t_scheme,sep="_")
        PrintDuration(ComputeDuration(t_init,getTime_()),
                      CallingFct)
    }
    
    list(Estim=Estim$Estim,duration=ComputeDuration(t_init,getTime_(),TRUE),
         method=method,tEstim=Estim$tEstim)
}

Compute2SGMMParametersEstim <- function(x,theta0,alphaReg,regularization,WeightingMatrix,
                                         t_scheme,pm,...){
    iter=0
    AllCurrentEstim <- ComputeCurrentEstim(t_scheme=t_scheme,theta0=theta0,x=x,alphaReg=alphaReg,
                                           regularization=regularization,WeightingMatrix="Id",
                                           pm=pm,...)
                                  
    theta1 <- (AllCurrentEstim$OptInfo)$par
    t <- AllCurrentEstim$t
    
    ProvidedWeightingMatrix <- ComputeWeightingMatrix(t,theta1,x,WeightingMatrix,pm,...)
    CurrentEstimOptInfo <- ComputeCurrentEstim(t_scheme=t_scheme,theta0=theta1,x=x,alphaReg=alphaReg,
                                               regularization=regularization,WeightingMatrix="provided",
                                               pm=pm,...,
                                               ProvidedWeightingMatrix=ProvidedWeightingMatrix)$OptInfo
    
    list(Estim=CurrentEstimOptInfo,tEstim=t)
}

ComputeITGMMParametersEstim <- function(x,theta0,alphaReg,regularization,WeightingMatrix,
                                        t_scheme,IterationControl,pm,...){
    iter=0
    Control <- checkIterationControl(IterationControl)
    AllCurrentEstim <-ComputeCurrentEstim(t_scheme=t_scheme,theta0=theta0,x=x,alphaReg=alphaReg,
                                          regularization=regularization,WeightingMatrix="Id",
                                          pm=pm,...)
    theta1 <-(AllCurrentEstim$OptInfo)$par
    t <- AllCurrentEstim$t
    PrevEstimParVal <- theta1
    RelativeErr=Control$RelativeErrMax+5
    
    while((iter < Control$NbIter) && (RelativeErr > Control$RelativeErrMax) ){
        ProvidedWeightingMatrix <- ComputeWeightingMatrix(t=t,theta=PrevEstimParVal,x=x,
                                                          WeightingMatrix=WeightingMatrix,
                                                          pm=pm,...)
        AllCurrentEstim <- ComputeCurrentEstim(t_scheme=t_scheme,
                                               theta0=PrevEstimParVal,x=x,
                                               alphaReg=alphaReg,
                                               regularization=regularization,
                                               WeightingMatrix="provided",
                                               pm=pm,...,
                                               ProvidedWeightingMatrix=ProvidedWeightingMatrix
                                               )
        CurrentEstimOptInfo <- AllCurrentEstim$OptInfo; t <- AllCurrentEstim$t
        CurrentEstimParVal <- CurrentEstimOptInfo$par

        if (Control$PrintIter) PrintIteration(CurrentEstimParVal,iter,Control$NbIter) 

        RelativeErr <- abs(CurrentEstimParVal-PrevEstimParVal)
        PrevEstimParVal <- CurrentEstimParVal
        iter=iter+1
    }
    list(Estim=CurrentEstimOptInfo,tEstim=t)
}

ComputeCueGMMParametersEstim <- function(x,theta0,alphaReg,regularization,WeightingMatrix,
                                         t_scheme,IterationControl,pm,...){
    iter=0
    Control <- checkIterationControl(IterationControl)
    PrevEstimParVal <- theta0
    RelativeErr=Control$RelativeErrMax+5
     
    while((iter < Control$NbIter) && (RelativeErr > Control$RelativeErrMax) ){
        AllCurrentEstim <- ComputeCurrentEstim(t_scheme=t_scheme,theta0=PrevEstimParVal,
                                               x=x,alphaReg=alphaReg,regularization=regularization,
                                               WeightingMatrix=WeightingMatrix,
                                               pm=pm,...)
        CurrentEstimOptInfo <- AllCurrentEstim$OptInfo; t <- AllCurrentEstim$t
        CurrentEstimParVal <- CurrentEstimOptInfo$par

        if (Control$PrintIter) PrintIteration(CurrentEstimParVal,iter,Control$NbIter)

        RelativeErr <- abs(CurrentEstimParVal-PrevEstimParVal)
        PrevEstimParVal <- CurrentEstimParVal
        iter=iter+1
    }

    list(Estim=CurrentEstimOptInfo,tEstim=t)
}

checkIterationControl <- function(IterationControl){
    NbIter <- ifelse(is.null(IterationControl$NbIter),
                     10,IterationControl$NbIter)
    PrintIter <- ifelse(is.null(IterationControl$PrintIter),
                        TRUE,IterationControl$PrintIter)
    RelativeErrMax <- ifelse(is.null(IterationControl$RelativeErrMax),
                             1e-3,IterationControl$RelativeErrMax)
    
    list(NbIter=NbIter,PrintIter=PrintIter,RelativeErrMax=RelativeErrMax)
}

ComputeCurrentEstim <- function(t_scheme,theta0,x,alphaReg,regularization,
                                WeightingMatrix,pm,...){

    t <- ComputeT(tScheme=t_scheme,theta=theta0,x=x,alphaReg=alphaReg,
                  regularization=regularization,
                  WeightingMatrix=WeightingMatrix,
                  pm=pm,...)
    
    optOutput <- nlminb(start=theta0,objective=ComputeObjective,
                        gradient = NULL, hessian = NULL,
                        t=t,x=x,alphaReg=alphaReg,
                        regularization=regularization,
                        WeightingMatrix=WeightingMatrix,
                        t_scheme=t_scheme,pm=pm,...,
                        lower = c(eps, -1+eps, 0+eps, -Inf),
                        upper = c(2-eps, 1-eps,  Inf,  Inf))
    
    list(OptInfo=optOutput,t=t)
}

ComputeObjective <- function(theta,t,x,alphaReg,regularization,
                             WeightingMatrix,t_scheme,pm,...){
    
    K <-ComputeWeightingMatrix(t=t,theta=theta,x=x,
                               WeightingMatrix=WeightingMatrix,
                               pm=pm,...)

    gbar <- sampleRealCFMoment(x=x,t=t,theta=theta,pm=pm)
    K1gbar <- ComputeInvKbyG(K=K,G=gbar,alphaReg=alphaReg,
                             regularization=regularization)
    
    obj <- crossprod(gbar,K1gbar)
    as.numeric(obj)
}

ComputeInvKbyG <- function(K,G,alphaReg,regularization){
    ComputeRegularized <- FALSE
    ## We compute Regularization if some eigenvalues are very small
    eigenAnalysis <- getSingularValueDecomposition(K)
    if (any(abs(eigenAnalysis$lambda) < alphaReg)) ComputeRegularized <- TRUE
    else{
        ## or if solve failed
        errStatut <- tryCatch(solve(K,G),error=function(e)e)
        err <- inherits(errStatut, "error")
        if (!err) K1G <- errStatut
        else ComputeRegularized <- TRUE
    }
    if (ComputeRegularized) K1G <- ComputeRegularizedK1G(K,G,alphaReg,regularization)

    K1G
}

ComputeRegularizedK1G <- function(K,G,alphaReg,regularization){
    RegularisedSol(Kn=K,alphaReg=alphaReg,r=G,regularization=regularization)
}

#This is actually the Inverse of the Variance: Called V in the report.
GMMasymptoticVarianceEstim <- function(...,t,theta,x,WeightingMatrix,
                                       alphaReg=0.1,regularization="Tikhonov",pm=0){
   
    K <- ComputeWeightingMatrix(t=t,theta=theta,x=x,WeightingMatrix=WeightingMatrix,pm=pm,...)
    B <- jacobianSampleRealCFMoment(t,theta,pm)

    fct <- function(G) ComputeInvKbyG(K=K,G=G,alphaReg=alphaReg,regularization=regularization)
    invKcrossB <- apply(X=B,MARGIN=2,FUN=fct)
   
    crossprod(B,invKcrossB)
}

## Note that we need:
## 1) Estimator convergence: n and n.alpha^(3/2) goes to inf, alpha goes to zero
## 2) Variance: n and n.alpha^(3) goes to inf, alpha goes to zero
test.GMMasymptoticVarianceEstim <- function(){
    WeightingMatrix="OptAsym"
    pm=0
    alphaReg=1
    theta <- c(1.5,0.5,1,0)
    t <- seq(0.1,3,length.out=100)
    n=10000; x=rstable(n,theta[1],theta[2],theta[3],theta[4],pm)

    V <- GMMasymptoticVarianceEstim(t=t,theta=theta,x=x,WeightingMatrix=WeightingMatrix,
                                    alphaReg=alphaReg)
    Var <- solve(V)/n
    list(V=V,V1=Var)
}

.asymptoticVarianceEstimGMM <- function(data,EstimObj,...){
    V <- solve(GMMasymptoticVarianceEstim(theta=EstimObj$Estim$par,
                                          t=EstimObj$tEstim,
                                          x=data,...))/length(data)                            
    NameParamsObjects(V)
}


getGMMmethodName <- function(algo,alphaReg,regularization,WeightingMatrix,
                             t_scheme,...){

    args <- list(...)
    if (!is.null(args$nb_t)) nt <- args$nb_t
    else if (!is.null(args$t_free)) nt <- length(args$t_free)
        
    paste(algo,
          paste("nb_t=",nt,sep=""),
          paste("alphaReg=",alphaReg,sep=""),
          paste("regularization=",regularization,sep=""),
          paste("WeightingMatrix=",WeightingMatrix,sep=""),
          paste("t_scheme=",t_scheme,sep=""),
          paste("OptimAlgo=","nlminb",sep=""),
          sep="_")
}
