## spacing: type of t spacing: Kout fro koutrelis type
##                           : UniformSpac: can use nb_t and/or nb_u
##                           : ArithSpac  : can use nb_t and/or nb_u                            
##                           : free: tpoints and upoints needs to be provided in the
## ...  args otherwise out type will be used

test.KoutParametersEstim<-function(){
    pm=0
    theta <- c(1.45,0.5,1.1,0.4)
    set.seed(1235);x <- rstable(500,theta[1],theta[2],theta[3],theta[4],pm=pm)
    theta0=theta-0.1
    spacing="Kout"#,"UniformSpac","ArithSpac","free")

    KoutParametersEstim(x=x,theta0=theta0,
                        spacing=spacing,pm=pm)
}

test.KoutObject <- function(){
    pm=0
    theta <- c(1.45,0.5,1.1,0.4)
    set.seed(1235);x <- rstable(500,theta[1],theta[2],theta[3],theta[4],pm=pm)
    Estim(EstimMethod="Kout",data=x,spacing="Kout",pm=pm)
}

getKoutMethodName <- function(spacing,...){
    paste("Koutrouvelis",
          paste("spacing=",spacing,sep=""),
          sep="_")
}

KoutParametersEstim<-function(x,theta0=NULL,
                              spacing=c("Kout","UniformSpac","ArithSpac","free"),
                              pm=0, tol=0.05,NbIter=10,
                              PrintTime=FALSE,...){

    spacing <- getSpacingType(match.arg(spacing),...)
    t_init<- getTime_()
    method <- getKoutMethodName(spacing)
    
    ## CHECK:Init initial Parameters with some methods
    if (is.null(theta0)) theta0 <- IGParametersEstim(x,pm)
    
    alpha <- theta0[1];beta <- theta0[2]
    gamma <- theta0[3];delta <- theta0[4]
    vals <- matrix(c(alpha,beta,gamma,delta),ncol=4)
    
    ## Shift data to fit the regression requirement
    if (pm==0) x <- x+beta*gamma*tan(pi2*alpha)

    ## Store old values
    alphaold <- alpha
    deltaold <- delta
    betaold <- beta
    deltaold <- delta
        
    for (iter in 1:NbIter)
        {
            s <- (x-delta)/gamma
            ParOld <- c(alpha,beta,gamma,delta)
            
            ## 1st regression: alpha and gamma
            reg1 <- EstimateAlphaGamma(RescaledData=s,
                                       spacing=spacing,
                                       oldPar=ParOld,
                                       iter=iter,...)
            alpha <- reg1$alpha ;gamma <- reg1$gamma
            s <- reg1$updatedData

            ## 2st regression : beta and delta
            reg2 <- EstimateBetaDelta(RescaledData=s,
                                       spacing=spacing,
                                       oldPar=c(alpha,beta,gamma,delta),
                                       iter=iter,...)
            beta <- reg2$beta ;delta <- reg2$delta
            s <- reg2$updatedData 
            
            vals <- rbind(vals,c(alpha,beta,gamma,delta))

            ## check for convergence
            diff=(alpha-ParOld[1])**2 + (delta-ParOld[4])**2
            if (is.na(beta)) beta <- betaold
            if (abs(beta)==1) par=c(alpha,betaold,gamma,deltaold)
            else par=c(alpha,beta,gamma,delta)
            if ((diff < tol)|| (abs(beta)==1)) break
        }
    
    if (PrintTime) PrintDuration(ComputeDuration(t_init,
                                                 t_final <- getTime_()
                                                 ),
                                 "KoutParametersEstim")
    res <- list(par=as.numeric(par),
                regObject1=reg1,regObject2=reg2,
                vals=vals)
    list(Estim=res,duration=ComputeDuration(t_init,getTime_(),TRUE),method=method)
}

.asymptoticVarianceEstimKout <- function(...){
    print("Kout Covarinace matrix not available")
    matrix(data=0,ncol=4,nrow=4)
}

getSpacingType <- function(spacing,...){
    args <- list(...)
    
    if (spacing=="free"){
        if (is.null(args$t_points) || is.null(args$u_points)) 
            spacing=="Kout" 
    }
    spacing
}

## first regression : Yk=log(-log{|phi(tk)|^2})=...

EstimateAlphaGamma <- function(RescaledData,spacing,oldPar,iter,...){
    s <- RescaledData
    alpha <- oldPar[1] ;beta <- oldPar[2];gamma <-oldPar[3]
    
    tp <- getTpoints(s,spacing,alpha,...)
    
    Phi <- .ecf(tp,s)
    y <- log(-log(Mod(Phi)^2))
    w <- log(Mod(tp))
    
    sigInit <-YCov(tp,length(s),alpha,beta,1) 
    sig <- checkCov(sigInit) 
    
    if(iter==1)reg1 <-lm(y~w)
    else reg1 <-lm.gls(y~w,W=sig,inverse=TRUE)

    m <- as.numeric((reg1[[1]])[1])
    alpha <- as.numeric((reg1[[1]])[2])
    invA <- 1/alpha
    gammahat <- as.numeric((0.5*exp(m))**invA)
    gamma=gamma*gammahat
    
    ## rescale and truncate
    s=s/gammahat
    alpha=max(alpha,0)
    alpha=min(alpha,2)
    gamma=max(gamma,0)

    list(alpha=alpha,gamma=gamma,
         obj=reg1,updatedData=s)
}

## 2nd regression : zi=arctan(Im(phi)/Re(phi))

EstimateBetaDelta <- function(RescaledData,spacing,oldPar,iter,...){
    s <- RescaledData
    alpha <- oldPar[1] ;beta <- oldPar[2];
    gamma <- oldPar[3];delta <- oldPar[4] 
    u <- getUpoints(s,spacing,alpha,...)
    
    
    z <- .get.arctan(s,u)
    Om <- (sign(u))*(abs(u)**alpha)
    sig2 <-ArcCov(u,length(s),alpha,beta,1)
    
    if (!(.checkPD(sig2))) {
        sig2 <- matrix(0,nrow=length(u),ncol=length(u))
        diag(sig2) <- diag(sig2+eps)
    }
    if(iter==1) reg2 <-lm(z~-1+u+Om)
    else reg2 <- lm.gls(z~-1+u+Om,W=sig2,inverse=TRUE)

    ## Update Parameters
    betaold <- beta
    beta <-  as.numeric((reg2[[1]])[2])/tan(pi2*alpha)
    Delta <- as.numeric((reg2[[1]])[1])
    deltaold <- delta;delta=Delta*gamma+delta #CHECK maybe delta <-> delta0
    s <- s-Delta
    
    ## rescale and truncate
    beta=min(beta,1);beta=max(beta,-1)
    list(beta=beta,delta=delta,
         obj=reg2,updatedData=s)
}


getTpoints <- function(data,spacing,alpha,...){
    if (is.null(list(...)$nb_t))
        nb_t <- getnbrTpoints(alpha,length(data))
    else nb_t <- list(...)$nb_t
    
    switch(spacing,
           "Kout"= {tt <- (pi/25)*seq(1:nb_t)},
           "UniformSpac"={
               An <- ComputeFirstRootRealeCF(x=data)
               tt <-seq(eps,An-eps,length.out=nb_t)
           },
           "ArithSpac"={
               An <- ComputeFirstRootRealeCF(x=data)
               tt <- (((1:nb_t)*(2:(nb_t+1)))/(nb_t*(nb_t+1)))*(An-eps)
           },
           "free"={tt <- list(...)$t_points}
           )
    tt
}

getUpoints <- function(data,spacing,alpha,...){
    if (is.null(list(...)$nb_u))
        nb_u <- getnbrLpoints(alpha,length(data))
    else nb_u <- list(...)$nb_u

    An <- ComputeFirstRootRealeCF(x=data)
    
    switch(spacing,
           "Kout"= {u <- seq(1:nb_u)*min(pi/50,An/nb_u)},
           "UniformSpac"={u <-seq(eps,An-eps,length.out=nb_u)},
           "ArithSpac"={
               u <- (((1:nb_u)*(2:(nb_u+1)))/(nb_u*(nb_u+1)))*(An-eps)
           },
           "free"={u <- list(...)$u_points}
           )
    u
}

YCov <- function (t ,N, alpha , beta,gam)
{
  # Compute covariance matrix of y = log (- log( phi(t) ) ), where phi(t) is 
  # ecf of alpha-stable random variables

    K = length(t)
    w = tan(alpha*pi/2)
    calpha = gam^alpha

    Tj = matrix(t ,nrow= 1 , ncol=K)
    Tk = matrix(t ,nrow= K , ncol=1)
    Tjalpha = abs(Tj)^alpha
    Tkalpha = abs(Tk)^alpha
    TkxTj = abs(Tk %*% Tj)
    TjpTk = matrix(Tj,ncol=K,nrow=K,byrow=T) + matrix(Tk,ncol=K,nrow=K)
    TjpTkalpha = abs(TjpTk)^alpha
    TjmTk = matrix(Tj,ncol=K,nrow=K,byrow=T) - matrix(Tk,ncol=K,nrow=K)
    TjmTkalpha = abs(TjmTk)^alpha

    A <- calpha*( matrix(Tjalpha,ncol=K,nrow=K,byrow=T) + matrix(Tkalpha,,ncol=K,nrow=K) - TjmTkalpha)
    B <- calpha * beta *
        (matrix(-Tjalpha * sign(Tj) * w,,ncol=K,nrow=K,byrow=T) 
        + matrix(Tkalpha * sign(Tk) * w,,ncol=K,nrow=K) 
        + TjmTkalpha * sign(TjmTk) * w) 
    D <- calpha * (matrix(Tjalpha,ncol=K,nrow=K,byrow=T) + matrix(Tkalpha,ncol=K,nrow=K) - TjpTkalpha)
    E <- calpha * beta *
        ( matrix(Tjalpha * sign(Tj) * w,ncol=K,nrow=K,byrow=T)    
        + matrix(Tkalpha * sign(Tk) * w,,ncol=K,nrow=K)
        - TjpTkalpha * sign(TjpTk) * w)
    
    res = (exp(A)* cos(B) + exp(D)*cos(E) - 2)/(2 * N * gam^(2*alpha) * TkxTj^alpha)
    res
 }

ArcCov <- function(t ,N, alpha , beta, gam)
{
  # Compute covariance matrix of z = Arctan(imag(phi(t))/real(phi(t)), 
  # where phi(t) is ecf of alpha-stable random variables.

    K = length(t)
    w = tan(alpha*pi/2)
    calpha = gam^alpha
    
    Tj = matrix(t ,nrow= 1 , ncol=K)
    Tk = matrix(t ,nrow= K , ncol=1)
    Tjalpha = abs(Tj)^alpha
    Tkalpha = abs(Tk)^alpha
    TkxTj = abs(Tk %*% Tj)
    TjpTk = matrix(Tj,ncol=K,nrow=K,byrow=T) + matrix(Tk,ncol=K,nrow=K)
    TjpTkalpha = abs(TjpTk)^alpha
    TjmTk = matrix(Tj,ncol=K,nrow=K,byrow=T) - matrix(Tk,ncol=K,nrow=K)
    TjmTkalpha = abs(TjmTk)^alpha

    B <- calpha * beta *
        (matrix(-Tjalpha * sign(Tj) * w, ncol=K,nrow=K,byrow=T)
        + matrix(Tkalpha * sign(Tk) * w ,ncol=K,nrow=K)
        + TjmTkalpha * sign(TjmTk) * w) 
    E <- calpha * beta *
        ( matrix(Tjalpha * sign(Tj) * w,ncol=K,nrow=K,byrow=T)    
        + matrix(Tkalpha * sign(Tk) * w,ncol=K,nrow=K)
        - TjpTkalpha * sign(TjpTk) * w)
    F <- calpha * (matrix(Tjalpha,ncol=K,nrow=K,byrow=T) + matrix(Tkalpha,ncol=K,nrow=K))
    G <- -calpha * TjmTkalpha
    H <- -calpha * TjpTkalpha
    
    res = exp(F) *(exp(G) * cos(B) - exp(H) * cos(E))/(2*N)
    res
  }
