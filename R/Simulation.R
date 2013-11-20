## CHECK: useful to get fct name
## For each couple of parameters in AlphaBetaMatrix:
## 1) if SaveEstim=TRUE :we save 1 file containing TruePar +Sample size value+seed +4 Parametes Value + Time + failure(if NaN or Na)
## 2) if StatSummary=TRUE : We produce a vector with :Sample size value+ alphaT,betaT + FctsToApply applied to the vector of estimation + Total nbr of failure + average time
## 3) Don't use both, prefer save file + function ComputeStatObjectFromFiles

Estim_Simulation <- function(AlphaBetaMatrix=abMat,SampleSizes=c(200,1600),MCparam=100,
                             Estimfct=c("ML","GMM","Cgmm","Kout"),HandleError=TRUE,
                             FctsToApply=StatFcts,saveOutput=TRUE,StatSummary=FALSE,
                             CheckMat=TRUE,tolFailCheck=tolFailure,
                             SeedOptions=NULL,...){
    
    SeedVector <- getSeedVector(MCparam,SeedOptions)
    Estimfct <- match.arg(Estimfct)
    nab <- nrow(AlphaBetaMatrix)
    lS <- length(SampleSizes)
    nRowOutput <- nab*lS
    indexStatOutput <- 1
    
    CheckPointValues <- readCheckPoint(AlphaBetaMatrix,Estimfct,nab,length(SampleSizes),MCparam,...)
    updatedCheckPointValues <- updateCheckPointValues(CheckPointValues,MCparam,lS,nab)
    if (updatedCheckPointValues$mc_start != 1 && StatSummary){
        print("Can't Compute Stat summary when the process doesn't the start from the beginning!!")
        StatSummary=FALSE
    }
    if (StatSummary){
        StatOutputLength <- length(FctsToApply)+ 5 # 5=Sample size+alphaT+betaT+Total failure+average time
        StatOutput <- list(alpha=matrix(data=NA,ncol=StatOutputLength,nrow=nRowOutput),
                           beta=matrix(data=NA,ncol=StatOutputLength,nrow=nRowOutput),
                           gamma=matrix(data=NA,ncol=StatOutputLength,nrow=nRowOutput),
                           delta=matrix(data=NA,ncol=StatOutputLength,nrow=nRowOutput))
    }
    
    for (ab in updatedCheckPointValues$ab_start:nab){
        alphaT <- AlphaBetaMatrix[ab,1]
        betaT <- AlphaBetaMatrix[ab,2]
        cat("---------------- a=",alphaT," *** b=",betaT," --------------- \n",sep="")

        if(saveOutput) initOutputFile(alphaT,betaT,MCparam,Estimfct,...)

        EstimOutput <- ComputeMCSimForAlphaBeta(alphaT=alphaT,betaT=betaT,MCparam=MCparam,
                                                SampleSizes=as.vector(SampleSizes),
                                                SeedVector=SeedVector,Estimfct=Estimfct,
                                                HandleError=HandleError,
                                                ab_current=ab,nab=nab,AlphaBetaMatrix,
                                                CheckPointValues=updatedCheckPointValues,
                                                saveOutput=saveOutput,...)
        if (StatSummary) {
            res <- ComputeStatOutput(EstimOutput=EstimOutput$outputMat,
                                     FctsToApply=FctsToApply,
                                     SampleSizes=SampleSizes,
                                     CheckMat=CheckMat,
                                     tolFailCheck=tolFailCheck,
                                     MCparam=MCparam,...)
            IndexSec <- seq(indexStatOutput,indexStatOutput+(lS-1),1)
            StatOutput$alpha[IndexSec,] <- res$alpha
            StatOutput$beta [IndexSec,] <- res$beta
            StatOutput$gamma[IndexSec,] <- res$gamma
            StatOutput$delta[IndexSec,] <- res$delta
            indexStatOutput <- indexStatOutput +lS
        }
    }
    deleteCheckPoint(AlphaBetaMatrix,Estimfct,nab,length(SampleSizes),MCparam,...)
    if (StatSummary) return(NameStatOutput(FctsToApply,StatOutput))
}

ComputeMCSimForAlphaBeta <- function(alphaT,betaT,MCparam,SampleSizes,
                                     SeedVector,Estimfct,HandleError,
                                     ab_current,nab,AlphaBetaMatrix,
                                     CheckPointValues,
                                     SaveOutput=TRUE,...){
    Ncol <- 10 # alphaT+betaT +Sample size +seed +4 Parametes Value + Time + failure
    nSS <- length(SampleSizes)
    Nrow <- nSS*MCparam
    Output <- matrix(data=NA,ncol=Ncol,nrow=Nrow)
    colnames(Output) <- c("alphaT","betaT","data size","seed","alphaE","betaE","gammaE","deltaE","failure","time")
    
    pm <- ifelse(is.null(args <- list(...)$pm),0,args)
    if (ab_current==CheckPointValues$ab_start){
        sample_start=CheckPointValues$sample_start
        mc_start=CheckPointValues$mc_start
    }
    else{
        sample_start=1
        mc_start=1
    }
    
    for (sample in sample_start:nSS){
        size <- SampleSizes[sample]
        if (sample != sample_start) mc_start=1
        for (mc in mc_start:MCparam){
            tIter <- getTime_();iter <- mc+(sample-1)*MCparam
            set.seed(seed <- SeedVector[mc])
            x <- rstable(n=size,alpha=alphaT,beta=betaT,
                         gamma=1,delta=0,pm=pm)

            Estim <- getEstimation(alphaT=alphaT,betaT=betaT,
                                   x=x,seed=seed,size=size,
                                   Ncol=Ncol,Estimfct=Estimfct,
                                   HandleError=HandleError,...)
            Output[iter,] <- Estim$outputMat ;file <- Estim$file
            ## update checkPoint and OutputFile
            writeCheckPoint(AlphaBetaMatrix,Estimfct,ab_current,nab,
                            sample,nSS,mc,MCparam,...)
            if (SaveOutput) updateOutputFile(alphaT,betaT,MCparam,Estim)   
            PrintEstimatedRemainingTime(iter,tIter,Nrow)
        }
    }
    list(outputMat=Output,file=file)
}

ComputeStatOutput <- function(EstimOutput,FctsToApply,SampleSizes,
                              CheckMat,tolFailCheck,MCparam,...){
    list(alpha=
         ComputeStatOutputPar(EstimOutput=EstimOutput,FctsToApply=FctsToApply,
                              par="alpha",SampleSizes=SampleSizes,
                              CheckMat=CheckMat,tolFailCheck=tolFailCheck,
                              MCparam=MCparam,...),
         beta=
         ComputeStatOutputPar(EstimOutput=EstimOutput,FctsToApply=FctsToApply,
                              par="beta",SampleSizes=SampleSizes,
                              CheckMat=CheckMat,tolFailCheck=tolFailCheck,
                              MCparam=MCparam,...),
         gamma=
         ComputeStatOutputPar(EstimOutput=EstimOutput,FctsToApply=FctsToApply,
                              par="gamma",SampleSizes=SampleSizes,
                              CheckMat=CheckMat,tolFailCheck=tolFailCheck,
                              MCparam=MCparam,...),
         delta=
         ComputeStatOutputPar(EstimOutput=EstimOutput,FctsToApply=FctsToApply,
                              par="delta",SampleSizes=SampleSizes,
                              CheckMat=CheckMat,tolFailCheck=tolFailCheck,
                              MCparam=MCparam,...)
         )
}

ComputeStatOutputPar <- function(EstimOutput,FctsToApply,par=c("alpha","beta","gamma","delta"),
                                 SampleSizes,CheckMat,tolFailCheck,MCparam,...){
    par <- match.arg(par)
    StatOutputLength <- length(FctsToApply)+ 5 # 5=Sample size+alphaT+betaT+Total failure+average time
    StatOutput=matrix(data=NA,ncol=StatOutputLength,nrow=length(SampleSizes))

    for (i in 1:length(SampleSizes)){
       n <- SampleSizes[i]
       EstimOutputOnseSampleSize <- getOneSampleEstim(SampleSize=n,par=par,
                                                      AllSampleEstimOutput=EstimOutput,
                                                      CheckMat=CheckMat,
                                                      tolFailCheck=tolFailCheck)
       StatOutput[i,] <- ComputeOutputOneSampleSize(i,FctsToApply,EstimOutputOnseSampleSize,par,MCparam,...)
    }
    StatOutput
}

getOneSampleEstim <- function(SampleSize,par,AllSampleEstimOutput,
                              CheckMat=TRUE,tolFailCheck=tolFailure){
    Allncol <- ncol(AllSampleEstimOutput)
    data <- AllSampleEstimOutput[AllSampleEstimOutput[,3]==SampleSize]
    AllMat <- matrix(data=data,ncol=Allncol)
    if (CheckMat) CheckedMat <- checkEstimMat(AllMat,tolFailCheck)
    else CheckedMat <- AllMat
        
    if(par=="alpha") mat <- matrix(data=CheckedMat[,c(-6,-7,-8)],ncol=Allncol-3)
    if(par=="beta")  mat <- matrix(data=CheckedMat[,c(-5,-7,-8)],ncol=Allncol-3)
    if(par=="gamma") mat <- matrix(data=CheckedMat[,c(-5,-6,-8)],ncol=Allncol-3)
    if(par=="delta") mat <- matrix(data=CheckedMat[,c(-5,-6,-7)],ncol=Allncol-3)
    mat     
}

checkEstimMat <- function(AllMat,tolFailCheck){
    NAVal <- rep(NA,4)
    res <- AllMat
    Tpar <- AllMat[1,1:2]
    
    fct <- function(v,Tpar) sum((v[5:6]-Tpar)**2) > tolFailCheck 
    check <- apply(AllMat,1,fct,Tpar=Tpar)
    res[check,5:8] <- NAVal
    res[check,ncol(AllMat)-1] <- 1
    
    res
}

ComputeOutputOneSampleSize <- function(row,FctsToApply,EstimOutput,par,MCparam,...){
    StatOutputLength <- length(FctsToApply)+ 5 # 5=Sample size+alphaT+betaT+Total failure+average time
    res <- matrix(data=NA,nrow=StatOutputLength,ncol=1)
    Tpar <- EstimOutput[1,1:2]
    n <- EstimOutput[1,3] # sample size
    
    if(row==1) res[1:2] <- Tpar
    res[3] <- n  # sample size
    res[4:(StatOutputLength-2)] <- ComputeFctsVals(FctsToApply,EstimOutput,par,...)
    #### \update Failure computation
    failCol <- EstimOutput[,ncol(EstimOutput)-1]
    res[StatOutputLength-1] <- sum(failCol) + (MCparam - length(failCol)) #failure
    res[StatOutputLength] <-mean(EstimOutput[,ncol(EstimOutput)]) # average time
    
    as.numeric(res)
}

ComputeFctsVals <- function(Fcts,EstimMat,par,...){
    pVal <- EstimMat[,5];pVal <- pVal[!is.na(pVal)]
    paramT <- c(EstimMat[1,1:2],1,0)
    n <- EstimMat[1,3]
    if(par=="alpha") index <- 1
    if(par=="beta")  index <- 2
    if(par=="gamma") index <- 3
    if(par=="delta") index <- 4
    
    res <- numeric(nf <- length(Fcts))
    for (i in 1:nf) res[i]=Fcts[[i]](p=pVal,paramT=paramT[index],n=n,...)
    
    res
}

## alphaT+betaT + Sample size +seed +4 Parametes Value + Time + failure
getEstimation <- function(alphaT,betaT,x,seed,size,Ncol,Estimfct,HandleError,...){
    output <- vector(length=Ncol)
    output[1:4] <- c(alphaT,betaT,size,seed)
    
    theta0 <- c(alphaT,betaT,1,0)-noise
    EstimRes <- Estim(EstimMethod=Estimfct,data=x,
                      theta0=theta0,ComputeCov=FALSE,
                      HandleError=HandleError,...)
    output[5:8] <- EstimRes@par
    output[9:10] <- c(EstimRes@failure,EstimRes@duration)
    list(outputMat=output,file=EstimRes@method)
}
