PrintIteration <- function(theta,iter,nbIterMax){
    print(paste("---------------iteration ", iter, "/",nbIterMax, "------------------------ (",
                theta[1],",",theta[2],",",theta[3],",",theta[4],")"))
    cat(" \n")
}

#==================================================================================================================
#----------------------------- Time -------------------------------------------------------------------------------
#==================================================================================================================
getTime_<- function() proc.time()[3]

ComputeDuration <- function(t_init,t_final,OneNumber=FALSE)
  {
    dt <- (t_final- t_init)
    dth <- floor(dt/3600)
    dtm <- floor((dt-dth*3600)/60)
    dts <- round((dt-dth*3600-dtm*60),0)
    if (OneNumber) t <- dt
    else t <- c(h=dth,min=dtm,sec=dts)
    t
  }

ConvertSecToVectorTime <- function(dt){
    dth <- floor(dt/3600)
    dtm <- floor((dt-dth*3600)/60)
    dts <- round((dt-dth*3600-dtm*60),0)

    c(h=dth,min=dtm,sec=dts)
}

ConvertVectorTimetoSec <- function(TimeVector){
    if (length(TimeVector)==3)
        dt <- TimeVector[1]*3600+TimeVector[2]*60+TimeVector[3]
    else dt <- as.numeric(TimeVector)
    dt
}

PrintDuration <- function(t,CallingFct="") {
    if (length(t)==1) tnew <- ConvertSecToVectorTime(t)
    else tnew <- t
    print(paste(CallingFct,":duration=",tnew[1]," h,",tnew[2]," min,",tnew[3]," sec.",paste=""))
}

## This function needs the starting time of the iteration procedure: a call to getTime()
PrintEstimatedRemainingTime <- function(ActualIter,ActualIterStartTime,TotalIterNbr){
    ti <- ConvertVectorTimetoSec(ActualIterStartTime)
    tf <- getTime_()

    AverageIterDur <- ComputeDuration(ti,tf,TRUE)
    RemainingDur <- (TotalIterNbr-ActualIter)*AverageIterDur
    dt <- ConvertSecToVectorTime(RemainingDur)
    cat("*** Iter ",ActualIter,"/",TotalIterNbr," *** Estimated Remaining Time: ",
        dt[1],"h",dt[2],"min",dt[3],"sec.",
        " *** \n",sep="")
}

#==================================================================================================================
#--------------------------------- Testing-------------------------------------------------------------------------
#==================================================================================================================
.almost_equal <- function(expected,tolExpect=1e-3){
    function(actual){
        if (is.list(actual) &&  is.list(expected) && length(actual)==length(expected)){
            diff <- 0
            for( i in 1:length(actual)){
                actualList <- actual[[i]]
                expectedList <- expected[[i]]
                diff <- diff+sum(abs(actualList-expectedList))
                answer <- as.logical(diff <= (length(actual)*tolExpect))
            }
        }
        else if (is.complex(actual) || is.complex(expected)){
            diff <- sum(Re(actual)-Re(expected))+sum(Im(actual)-Im(expected))
            answer <- as.logical(abs(diff) <= (2*tolExpect))
        }
        else{
            answer <- as.logical((diff <- sum(abs(actual-expected))) <= tolExpect)
        }


        type <- identical(answer, TRUE)
        ## 2016-07-26 TODO: make this unconditiona in the near future.
        if(packageVersion("testthat") >= "1.0.0")
            ## should be in c("success", "failure", "error", "skip", "warning")
            type <- if(type) "success"
                    else "failure"

        expectation(type,
                    paste("absolute error= ",diff ," > ", tolExpect))
    }
}

expect_almost_equal <- function(x,y,tolExpect=1e-3){
    expect_that(x,.almost_equal(y,tolExpect))}

#===================================================================================================================
#----------------------------------- Computation -------------------------------------------------------------------
#===================================================================================================================

## Compute the cross sum of 2 vectors
## Return a matrix nrow(X) x nrow(Y)
## Res[i,j]=x(i)+y(j)

crossSum <- function(X,Y){
    fct <- function(x) X+x
    res <- sapply(Y,fct)
    res
}

## the choice of 3+ seet Initial value is arbitrary

getSeedVector <- function(Outputsize,SeedOptions){
    set.seed(345)
    if (is.null(SeedOptions))
        vec <- as.vector(sample.int(n=3*Outputsize,size=Outputsize))
    else {
        MCtot <- SeedOptions$MCtot
        seedStart <- SeedOptions$seedStart
        seedEnd <- seedStart+Outputsize
        vec <- as.vector(sample.int(n=3*MCtot,size=MCtot))[seedStart:seedEnd]
    }
    vec
}

#===================================================================================================================
#---------------------------------- Formatting ---------------------------------------------------------------------
#===================================================================================================================

NameParamsObjects <- function(mat){
    parNames <- c("alpha","beta","gamma","delta")
    minMaxCol <- c("min","max")

    if (length(mat)==4){
        names(mat) <- parNames
    }
    else if (is.matrix(mat) && nrow(mat)==4){
        rownames(mat) <- parNames
        if (ncol(mat)==2) colnames(mat) <-minMaxCol
        else if (ncol(mat)==4) colnames(mat) <- parNames
    }
    mat
}

## if the matrix if bad conditioned, we set to 0
## the off diagonal coefficient.

checkCov <- function(sig)
    {
        res <- matrix(0,ncol=ncol(sig),nrow=nrow(sig))
        crit_cond <- 9.9e-15 # critical value of the condition number
        condSig <- -Inf

        tr <- tryCatch(rcond(sig),
                       error=function(e)e)
        err <- inherits(tr, "error")
        if ( (!err) &&  ( !is.nan(tr) ) && ( !is.na(tr) ) ) {condSig <- tr}

        if (condSig < crit_cond){
            diag(res) <- diag(sig)
        }
        else res <-sig
        res
    }


## Check if matrix W is definite Positive by computing the eigenvalues

.checkPD <- function(W){
  eW <- eigen(W, TRUE)
  d <- eW$values
  res <- ifelse (any(d <= 0),FALSE,TRUE)
  res
}

## Compute the arctan:
## smooth the gn in order to make sure all values
## accounts for non principal branches

.get.arctan <- function(x,u)
  {
    PhiU <- .ecf(u,x)
    g <- atan2(Im(PhiU),Re(PhiU))
    res <- smooth_gu(g)
    res
  }

## function to take into account the possible nonprincipal
## branches of arctan function

smooth_gu <- function(g)
    {
        m <- min(g);M <- max(g)
        if ((sup <- (min(g)+pi)) > M) res <- g # fct ctn

        else {
            res <- numeric(length(g))
            for (i in 1:length(g))
                {
                    res[i] <- g[i]
                    while (res[i] > sup){
                        res[i] <- res[i]-pi
                    }
                }
        }
        res
    }
