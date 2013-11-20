phinR <- function(t,x) mean(cos(t*x))

ComputeFirstRootRealeCF <- function(x,...,tol=1e-3,maxIter=100,
                                    lowerBand=1e-4,upperBand=30){
    
    WelshSol <- WelshFirstRootRealeCF(x,tol,maxIter)
    if (WelshSol$phinR < tol) return(WelshSol$t)
    else return(numFirstRootRealeCF(x,tol,lowerBand,upperBand,...)$t)
}

## Based on the procedure suggested by A.H WELSH(1986) in
## Implemmenting empirical characteristic function procedures

WelshFirstRootRealeCF <- function(x,tol=1e-3,maxIter=100){
    A=0;iter=0
    m=mean(abs(x))
    val=phinR(A,x)
    
    while ((abs(val) > tol) && (iter< maxIter)){
        A=A+val/m
        val=phinR(A,x)
        iter=iter+1
    }
    list(t=A,phinR=val)
}

graphFirstRootRealeCF <- function(x,tol=1e-3,lowerBand=1e-4,upperBand=30){
    t_seq<- seq(lowerBand,upperBand,tol)
    phiVal <- sapply(t_seq,phinR,x=x)
    
    t <- t_seq[abs(phiVal)< tol][1]
    list(t=t, phinR=phinR(t,x))
}

numFirstRootRealeCF <- function(x,tol=1e-3,lowerBand=1e-4,upperBand=30,...){
    t_init<-graphFirstRootRealeCF(x,tol=tol,
                                  lowerBand=lowerBand,
                                  upperBand=upperBand)$t
    
    if (is.na(t_init)) t_init <- upperBand
    objectiveFct <- function(t) abs(phinR(t,x))
    
    optInfo <- nlminb(start=t_init,objective=objectiveFct,
                      lower=lowerBand,
                      upper=upperBand)
    
    list(t=as.numeric(optInfo$par),phinR=optInfo$objective)
}

test.ComputeComputeFirstRootRealeCF <- function(){
    test.WelshFirstRootRealeCF()
    test.graphFirstRootRealeCF()
    test.numFirstRootRealeCF()
}

test.numFirstRootRealeCF <- function(){
    set.seed(345); x <- rstable(500,1.5,0.5)
    tEstim <- numFirstRootRealeCF(x)$t
    tRef <-  2.305364
    expect_almost_equal(tEstim,tRef)
}

test.graphFirstRootRealeCF <- function(){
    set.seed(345); x <- rstable(500,1.5,0.5)
    tEstim <- graphFirstRootRealeCF(x)$t
    tRef <- 2.3031
    expect_almost_equal(tEstim,tRef)
}

test.WelshFirstRootRealeCF <- function(){
    set.seed(345); x <- rstable(500,1.5,0.5)
    tEstim <- WelshFirstRootRealeCF(x)$t
    tRef <- 2.302698
    expect_almost_equal(tEstim,tRef)
}
