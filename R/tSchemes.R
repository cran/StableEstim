ComputeT <- function(tScheme,...){
    args <- list(...)
    
    if(tScheme=="free") {
        checkFreeArgs(args)
        t <- args$t_free
    }
    else if((tScheme=="equally") || (tScheme=="NonOptAr")){
        checkMainArgs(args)
        t <- ComputeEquallySpacedPoints(tScheme=tScheme,...)
    }
    else if((tScheme=="uniformOpt") || (tScheme=="ArithOpt")
       || (tScheme=="VarOpt")){
        nb_t <- checkMainArgs(args)
        t <- ComputeOptimisedPoints(tScheme=tScheme,...)                                   
    }
    t
}

checkFreeArgs <- function(args){
    if (is.null(args$t_free))
        stop("You need to provide t when you choose free as tScheme")
}

checkMainArgs <- function(args){
    if(is.null(args$x)) stop("you didn't provide x to compute t")
    if (is.null(args$nb_t))
        warning("You didn't provide a number of points to compute t, we will use 40")
}

ComputeEquallySpacedPoints <- function(...,tScheme,x,nb_t,Constrained=TRUE){
    Bands <- Compute_tBands(x=x,Constrained=Constrained,...)
    if (tScheme=="equally"){ 
        t <- seq(Bands$lower,Bands$upper,length.out=nb_t)
    }
    else if (tScheme=="NonOptAr"){
       t <- (((1:nb_t)*(2:(nb_t+1)))/(nb_t*(nb_t+1)))*Bands$upper
    }
    t
}

Compute_tBands <- function(x,Constrained,...){
    args <- list(...)
    if (!Constrained){
        lower <- ifelse(is.null(args$min_t),eps,args$min_t)
        upper <- ifelse(is.null(args$max_t),ComputeFirstRootRealeCF(x=x,...)-eps,
                        args$max_t)
    }
    else{
        An <- ComputeFirstRootRealeCF(x=x,...)
        lower=eps; upper=An-eps
    }
    list(lower=lower,upper=upper)
}

## Constrained =TRUE if you want to find t in [eps,An]
ComputeOptimisedPoints <- function(...,tScheme,nb_t=40,Constrained=TRUE,
                                   FastOptim=TRUE){
    An <- ComputeFirstRootRealeCF(...)
    t0 <- seq(eps,An-eps,length.out=nb_t)
    
    if((tScheme=="uniformOpt") || (tScheme=="ArithOpt")){
        t <- ComputeApproxOptSpacing(...,tScheme=tScheme,nb_t=nb_t,
                                     Constrained=Constrained,An=An)
    }
    else if (tScheme=="VarOpt"){
        if (FastOptim)
            t <- optim(par=t0,fn=ObjectiveFctToMinIn_t,
                       gr=NULL,...,
                       method = "Nelder-Mead")$par
        else
            t <- nlminb(start=t0,objective=ObjectiveFctToMinIn_t,...)
    }
    t
}

ComputeApproxOptSpacing <- function(...,tScheme,nb_t,Constrained,An){
    tau0 <- ComputeTau0(An=An,tScheme=tScheme,nb_t=nb_t)
    if (Constrained) {
        Bands <- Compute_tauBands(tScheme,nb_t,An)
        tauInfo <- nlminb(start=tau0,objective=ObjectiveFctToMinIn_tau,
                          gradient = NULL, hessian = NULL,
                          tScheme=tScheme,nb_t=nb_t,...,
                          lower=Bands$lower,upper=Bands$upper)
    }
    else{
        tauInfo <- nlminb(start=tau0,objective=ObjectiveFctToMinIn_tau,
                          gradient = NULL, hessian = NULL,
                          tScheme=tScheme,nb_t=nb_t,...)
    }
    tau <- tauInfo$par
    if(tScheme=="uniformOpt"){t <- (1:nb_t)*tau}
    else if (tScheme=="ArithOpt") {t <- ((1:nb_t)*(2:(nb_t+1)))*tau}
    t    
}

ComputeTau0 <- function(An,tScheme,nb_t){
    if(tScheme=="uniformOpt")
        tau0 <- An/(2*nb_t)
    else if (tScheme=="ArithOpt") 
        tau0 <- (An-eps)/(nb_t*(2*nb_t+1))
}

Compute_tauBands <- function(tScheme,nb_t,An){
    if(tScheme=="uniformOpt"){
        lower=eps
        upper=An/nb_t
    }
    else if (tScheme=="ArithOpt"){ 
        lower=min(eps/2,An/(nb_t*(nb_t**2+1)))
        upper=(An-eps)/(nb_t*(nb_t+1))
    }
    Bands <- list(lower=lower,upper=upper)
}

## det(A^-1)=1/det(A) ;
## abs is necessary to avoid any sign problem ( we want a det close to zero)
InvDet <- function(Mat) 1/abs(det(Mat))

ObjectiveFctToMinIn_tau <- function(tau,...,tScheme,nb_t,theta,x,WeightingMatrix,
                                  alphaReg=0.01,regularization="Tikhonov",
                                  pm=0,fctToApply=InvDet){

    if(tScheme=="uniformOpt"){t <- (1:nb_t)*tau}
    else if (tScheme=="ArithOpt") {t <- ((1:nb_t)*(2:(nb_t+1)))*tau} 

    ObjectiveFctToMinIn_t(t,...,theta=theta,x=x,WeightingMatrix=WeightingMatrix,
                          alphaReg=alphaReg,regularization=regularization,
                          pm=pm,fctToApply=fctToApply)
}

## if V is the variance of the GMM estimator, then V=W^-1
## where W=B'(K^-1)B.
## det(V)=1/det(W)
ObjectiveFctToMinIn_t <- function(t,...,theta,x,WeightingMatrix,
                                  alphaReg=0.01,regularization="Tikhonov",
                                  pm=0,fctToApply=InvDet){
    
    W <- GMMasymptoticVarianceEstim(...,t=t,theta=theta,x=x,WeightingMatrix=WeightingMatrix,
                                    alphaReg=alphaReg,regularization=regularization,pm=pm)
    fctToApply(W)
}
