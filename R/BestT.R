## ... is to be passed to ComputeWeightingMatrix or optim
## to need it here

ComputeBest_t <- function(AlphaBetaMatrix=abMat,nb_ts=seq(10,100,10),alphaReg=0.001,
                          FastOptim=TRUE,...){
    output <- list()
    nab <- nrow(AlphaBetaMatrix)
        
    for (ab in 1:nab){
        theta <- c(AlphaBetaMatrix[ab,1],AlphaBetaMatrix[ab,2],1,0)
        cat("---------------- a=",theta[1]," *** b=",theta[2]," --------------- \n",sep="")
        OnethetaObj <- new(Class="Best_t",theta=theta,nbt=nb_ts)
        for (it in 1:length(nb_ts)){
            nt <- nb_ts[it]
            An <- 4
            t0 <- seq(eps,An-eps,length.out=nt)
            cat("---------------- nt=",nt," -------------------- \n",sep="")
            if (FastOptim){
                tAll <- optim(par=t0,fn=ObjectiveFctToMinIn_t,
                              theta=theta,x=NULL,
                              WeightingMatrix="OptAsym",alphaReg=alphaReg,...,
                              method="Nelder-Mead")
                OnethetaObj@detVal[it] <- as.numeric(tAll$value)
            }
            else{
                tAll <- nlminb(start=t0,objective=ObjectiveFctToMinIn_t,
                               theta=theta,x=NULL,
                               WeightingMatrix="OptAsym",alphaReg=alphaReg,...)
                OnethetaObj@detVal[it] <- as.numeric(tAll$objective)
            }
            OnethetaObj@tvec[[it]] <- as.numeric(tAll$par)
            OnethetaObj@convcode[it] <- as.numeric(tAll$convergence)
        }
        output[[ab]] <- OnethetaObj 
    }
    output
}

ComputeBest_tau <- function(AlphaBetaMatrix=abMat,nb_ts=seq(10,100,10),
                            tScheme=c("uniformOpt","ArithOpt"),
                            Constrained=TRUE,alphaReg=0.001,...){
    output <- list()
    nab <- nrow(AlphaBetaMatrix)

    for (ab in 1:nab){
        theta <- c(AlphaBetaMatrix[ab,1],AlphaBetaMatrix[ab,2],1,0)
        cat("---------------- a=",theta[1]," *** b=",theta[2]," --------------- \n",sep="")
        OnethetaObj <- new(Class="Best_t",theta=theta,nbt=nb_ts)
        for (it in 1:length(nb_ts)){
            nt <- nb_ts[it]
            An <- 4
            tau0 <- ComputeTau0(An=An,tScheme=tScheme,nb_t=nt)
            if (Constrained)
                Bands <- Compute_tauBands(tScheme,nt,An)
            else Bands <- list(lower=-Inf,upper=Inf)
            
            cat("---------------- nt=",nt," -------------------- \n",sep="")
            tau_All <- nlminb(start=tau0,objective=ObjectiveFctToMinIn_tau,
                              theta=theta,x=NULL,
                              tScheme=tScheme,nb_t=nt,
                              WeightingMatrix="OptAsym",
                              alphaReg=alphaReg,...,
                              lower=Bands$lower,upper=Bands$upper)
            
            OnethetaObj@tvec[[it]] <- .getTfromTau(as.numeric(tau_All$par),
                                                  tScheme,
                                                  nt)
            OnethetaObj@detVal[it] <- as.numeric(tau_All$objective)
            OnethetaObj@convcode[it] <- as.numeric(tau_All$convergence)
        }
        output[[ab]] <- OnethetaObj
    }
    output
}

.getTfromTau <- function(tau,tScheme=c("uniformOpt","ArithOpt"),nt){
    if (tScheme=="uniformOpt") t <- (1:nt)*tau
    else if (tScheme=="ArithOpt") t <- ((1:nt)*(2:(nt+1)))*tau
    t
}
