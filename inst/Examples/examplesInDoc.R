example.GMM <- function(){
    ## General data
    theta <- c(1.5,0.5,1,0)
    pm <- 0
    set.seed(345);
    x <- rstable(100,theta[1],theta[2],theta[3],theta[4],pm)
    ##---------------- 2S free ----------------
    ## method specific arguments
    regularization="cut-off"
    WeightingMatrix="OptAsym"
    alphaReg=0.005
    
    ## If you are just interested by the value
    ## of the 4 estimated parameters
    t_scheme="free"
    t_seq=seq(0.1,2,length.out=12)
    algo = "2SGMM"
    
    suppressWarnings(GMMParametersEstim(x=x,
                                        algo=algo,alphaReg=alphaReg,
                                        regularization=regularization,
                                        WeightingMatrix=WeightingMatrix,
                                        t_scheme=t_scheme,
                                        pm=pm,PrintTime=TRUE,t_free=t_seq))   
}

example.ML <- function(){
    theta <- c(1.5,0.4,1,0)
    
    pm <- 0
    ## 50 points does not give accurate estimation
    ## but it makes estimation fast for installation purposes
    ## use at least 200 points to get decent results.
    set.seed(1333);x <- rstable(50,theta[1],theta[2],theta[3],theta[4],pm)

    MLParametersEstim(x=x,pm=pm,PrintTime=TRUE)

}

example.Cgmm <- function(){
    ## general inputs
    theta <- c(1.45,0.55,1,0)
    pm <- 0
    set.seed(2345);x <- rstable(200,theta[1],theta[2],theta[3],theta[4],pm)

    ## GMM specific params
    alphaReg=0.01
    subdivisions=20
    randomIntegrationLaw="unif"
    IntegrationMethod="Simpson"
    IterationControl <- list(NbIter=3,PrintIter=FALSE,RelativeErrMax=1e-2)
    
    ## Estimation
    CgmmParametersEstim(x=x,type="IT",alphaReg=alphaReg,
                        subdivisions=subdivisions,
                        IntegrationMethod=IntegrationMethod,
                        randomIntegrationLaw=randomIntegrationLaw,
                        s_min=0,s_max=1,
                        IterationControl=IterationControl,
                        pm=pm,PrintTime=TRUE)
}

example.Estim <- function(){
    theta <- c(1.45,0.55,1,0)
    pm <- 0
    set.seed(2345)
    x <- rstable(200,theta[1],theta[2],theta[3],theta[4],pm)

    ## Cgmm specific inputs 
    alphaReg=0.01
    theta0 <- c(1.35,0.45,1,0)
    subdivisions=20
    
    randomIntegrationLaw="unif"
    IntegrationMethod="Simpson"
    
    ## Estim procedure
    Estim(EstimMethod="Cgmm",data=x,
          ComputeCov=TRUE,HandleError=FALSE,
          type="2S",alphaReg=alphaReg,
          subdivisions=subdivisions,
          IntegrationMethod=IntegrationMethod,
          randomIntegrationLaw=randomIntegrationLaw,
          s_min=0,s_max=1,theta0=theta0,pm=pm,
          PrintTime=TRUE,level=0.92)

}

example.RegSol <- function(){
    ## Adapted from R examples for Solve 
    ## We compare the result of the regularized sol to the expected solution
    
    hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+")}
    
    K_h8 <- hilbert(8);
    r8 <- 1:8
    
    alphaReg_robust<- 1e-4
    Sa8_robust <- RegularisedSol(K_h8,alphaReg_robust,r8,"LF")
    
    alphaReg_accurate<- 1e-10
    Sa8_accurate <- RegularisedSol(K_h8,alphaReg_accurate,r8,"LF")
    
    ## when pre multiplied by K_h8 ,the expected solution is 1:8
    ## User can check the influence of the choice of alphaReg

    list(Sa8_robust,Sa8_accurate)
}

example.Integrate <- function(){
    ## Define the integrand
    f_fct <- function(s,x){sapply(X=x,
                                  FUN=sampleComplexCFMoment,
                                  t=s,theta=theta)
                       }
    f_bar_fct <- function(s,x){Conj(f_fct(s,x))}
    
    ## Function specific arguments
    theta <- c(1.5,0.5,1,0)
    set.seed(345);X=rstable(3,1.5,0.5,1,0)
    s_min=0;s_max=2
    numberIntegrationPoints=10
    randomIntegrationLaw="norm"
    
    IntegrateRandomVectorsProduct(f_fct,X,f_bar_fct,X,s_min,s_max,
                                  numberIntegrationPoints,
                                  "Simpson",randomIntegrationLaw)
}

example.cf <- function(){
    ## define the parameters
    nt <- 10
    t <- seq(0.1,3,length.out=nt)
    theta <- c(1.5,0.5,1,0)
    pm <- 0
    
    ## Compute the characteristic function
    ComplexCF(t=t,theta=theta,pm=pm)
}

example.jack <- function(){
    ## define the parameters
    nt <- 10
    t <- seq(0.1,3,length.out=nt)
    theta <- c(1.5,0.5,1,0)
    pm <- 0
    
    ## Compute the jacobian of the characteristic function
    jacobianComplexCF(t=t,theta=theta,pm=pm)
}

example.Kout <- function(){
    pm=0
    theta <- c(1.45,0.5,1.1,0.4)
    set.seed(1235);x <- rstable(500,theta[1],theta[2],theta[3],theta[4],pm=pm)
    theta0=theta-0.1
    spacing="Kout"

    KoutParametersEstim(x=x,theta0=theta0,
                        spacing=spacing,pm=pm)


}

example.IG <- function(){
     x <- rstable(200,1.2,0.5,1,0,pm=0)
     IGParametersEstim(x,pm=0)
}
