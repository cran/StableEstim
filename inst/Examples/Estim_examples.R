## -----------------------------------  estimation parameters -------------------------
theta <<- c(1.45,0.55,1,0)
pm <<- 0
set.seed(2345);x <<- rstable(1000,theta[1],theta[2],theta[3],theta[4],pm)
theta0 <<- c(1.35,0.45,1,0)
##------------------------------------ Kout ------------------------------------------
example.Kout <- function(){
    objKout <- Estim(EstimMethod="Kout",data=x,
                     ComputeCov=FALSE,HandleError=FALSE,
                     spacing="Kout")
    
    objKoutUniform <- Estim(EstimMethod="Kout",data=x,
                            ComputeCov=FALSE,HandleError=FALSE,
                            spacing="ArithSpac",nb_t=23)
    
    objKoutArith <- Estim(EstimMethod="Kout",data=x,
                          ComputeCov=FALSE,HandleError=FALSE,
                          spacing="UniformSpac",nb_u=13)
    
    objKoutFree <- Estim(EstimMethod="Kout",data=x,
                         ComputeCov=FALSE,HandleError=FALSE,
                         spacing="free",t_points=seq(0.2,2.5,length.out=12),
                         u_points=seq(0.1,1.85,length.out=10))
    
    list(kout=objKout,
         Uniform=objKoutUniform,
         Arith=objKoutArith,
         Free=objKoutFree
         )
}
##------------------------------------ Cgmm ------------------------------------------
example.Cgmm <- function(index=1){
    subdivisions=50
    alphaReg=0.01
    randomIntegrationLaw="unif"
    IntegrationMethod="Simpson"

    if (index==1){ ## Takes around 12 min
        obj2SCgmm <- Estim(EstimMethod="Cgmm",data=x,
                           ComputeCov=TRUE,HandleError=FALSE,
                           type="2S",alphaReg=alphaReg,
                           subdivisions=subdivisions,
                           IntegrationMethod=IntegrationMethod,
                           randomIntegrationLaw=randomIntegrationLaw,s_min=0,s_max=1,
                           theta0=NULL,pm=pm,PrintTime=FALSE,level=0.95)
        objITCgmm <- NULL; objCueCgmm <- NULL
    }
    else { 
        IterationControl <- list(NbIter=3,PrintIter=FALSE,RelativeErrMax=1e-2)
        randomIntegrationLaw="norm"
        IntegrationMethod="Uniform"
        
        if (index==2){ # takes 6 min
            objITCgmm <- Estim(EstimMethod="Cgmm",data=x,
                               ComputeCov=TRUE,HandleError=FALSE,
                               type="IT",alphaReg=alphaReg,
                               subdivisions=subdivisions,
                               IntegrationMethod=IntegrationMethod,
                               randomIntegrationLaw=randomIntegrationLaw,s_min=0,s_max=1,
                               theta0=theta0,IterationControl=IterationControl,pm=pm,
                               PrintTime=FALSE,level=0.92,mu=0.5,sigma=2) # mu and sigma: params of the normal law
            obj2SCgmm <- NULL ;objCueCgmm <- NULL
        }
        else { #takes 17 min
            objCueCgmm <- Estim(EstimMethod="Cgmm",data=x,
                                ComputeCov=TRUE,HandleError=FALSE,
                                type="Cue",alphaReg=alphaReg,
                                subdivisions=subdivisions,
                                IntegrationMethod=IntegrationMethod,
                                randomIntegrationLaw=randomIntegrationLaw,s_min=0,s_max=1,
                                theta0=theta0,IterationControl=IterationControl,pm=pm,
                                PrintTime=FALSE,level=0.92)
            obj2SCgmm <- NULL ;objITCgmm <- NULL
        }
    }
    
    list(twoStep=obj2SCgmm,
         IT=objITCgmm,
         Cue=objCueCgmm)
}

##------------------------------------ ML ------------------------------------------
example.ML <- function(){# takes 4min
    Estim(EstimMethod="ML",data=x,theta0=NULL,ComputeCov=TRUE,HandleError=FALSE,
          level=0.90,subdivisions=50,pm=0)
}

##----------------------------------- GMM ------------------------------------------
example.GMM.Classic <- function(){ #takes less than 1 min
        ## free
    t_scheme="free"
    algo <- "2SGMM"
    alphaReg=0.005
    regularization="Tikhonov"
    WeightingMatrix="DataVar"
    t_seq=seq(0.1,2,length.out=42)
    
    objfree <- Estim(EstimMethod="GMM",data=x,theta0=theta0,
                     ComputeCov=TRUE,HandleError=FALSE,
                     algo=algo,alphaReg=alphaReg,
                     regularization=regularization,
                     WeightingMatrix=WeightingMatrix,
                     t_scheme=t_scheme,
                     IterationControl=IterationControl,
                     PrintTime=TRUE,
                     pm=pm,t_free=t_seq,level=0.92)
    
    ## equally
    IterationControl <- list(NbIter=8,PrintIter=FALSE,RelativeErrMax=1e-2)
    t_scheme="equally"
    algo="ITGMM"
    regularization="Tikhonov"
    WeightingMatrix="DataVar"
    
    objequally <- Estim(EstimMethod="GMM",data=x,theta0=theta0,
                        ComputeCov=TRUE,HandleError=FALSE,
                        algo=algo,alphaReg=alphaReg,
                        regularization=regularization,
                        WeightingMatrix=WeightingMatrix,
                        t_scheme=t_scheme,
                        IterationControl=IterationControl,
                        pm=pm,
                        PrintTime=TRUE,nb_t=42,Constrained=TRUE,
                        tol=1e-3,maxIter=100,lowerBand=0.01,
                        upperBand=10,min_t=0.1,max_t=3)
    
    ## NonOptAr
    t_scheme="NonOptAr"
    algo="CueGMM"    
    regularization="LF"
    WeightingMatrix="OptAsym"
    
    objNonOptAr <- Estim(EstimMethod="GMM",data=x,theta0=theta0,
                         ComputeCov=TRUE,HandleError=FALSE,
                         algo=algo,alphaReg=alphaReg,
                         regularization=regularization,
                         WeightingMatrix=WeightingMatrix,
                         t_scheme=t_scheme,
                         IterationControl=IterationControl,
                         pm=pm,
                         PrintTime=TRUE,nb_t=42,Constrained=FALSE,
                         tol=1e-3,maxIter=100,lowerBand=0.01,
                         upperBand=10,min_t=0.1)

    list(NonOptAr=objNonOptAr,equally=objequally,free=objfree)
}

example.GMM.opt <- function(index=1){
    theta <- c(1.13,0.23,1,0)
    theta0 <- NULL
    pm <- 0
    set.seed(345);x <- rstable(500,theta[1],theta[2],theta[3],theta[4],pm)
    alphaReg=0.05
    IterationControl <- list(NbIter=8,PrintIter=FALSE,RelativeErrMax=1e-2)
    
    if(index==1){##takes  2s
        ## uniformOpt
        t_scheme="uniformOpt"
        algo="CueGMM"
        regularization="cut-off"
        WeightingMatrix="Id"
    
        obj <- Estim(EstimMethod="GMM",data=x,
                     ComputeCov=TRUE,HandleError=FALSE,
                     algo=algo,alphaReg=alphaReg,
                     regularization=regularization,
                     WeightingMatrix=WeightingMatrix,
                     t_scheme=t_scheme,
                     IterationControl=IterationControl,
                     pm=pm,
                     PrintTime=TRUE,nb_t=42,Constrained=TRUE,
                     tol=1e-3,maxIter=100,lowerBand=0.01,
                     upperBand=10)
    }
    else if (index==2){ ## takes 4s
        ## ArithOpt
        t_scheme="ArithOpt"
        algo="ITGMM"
        regularization="LF"
        WeightingMatrix="DataVar"
        nb <- 42
        
        obj <- Estim(EstimMethod="GMM",data=x,theta0=theta0,
                     ComputeCov=TRUE,HandleError=FALSE,
                     algo=algo,alphaReg=alphaReg,
                     regularization=regularization,
                     WeightingMatrix=WeightingMatrix,
                     t_scheme=t_scheme,
                     IterationControl=IterationControl,
                     pm=pm,
                     PrintTime=TRUE,nb_t=42,Constrained=TRUE,
                     tol=1e-3,maxIter=100,lowerBand=0.01,
                     upperBand=10)
    }
    else if (index==3){ #takes 22s
        ## VarOpt
        t_scheme="VarOpt"
        algo="2SGMM"    
        regularization="Tikhonov"
        WeightingMatrix="OptAsym"
        nb <- 42
        
        obj <- Estim(EstimMethod="GMM",data=x,theta0=theta0,
                     ComputeCov=TRUE,HandleError=FALSE,
                     algo=algo,alphaReg=alphaReg,
                     regularization=regularization,
                     WeightingMatrix=WeightingMatrix,
                     t_scheme=t_scheme,
                     IterationControl=IterationControl,
                     pm=pm,
                     PrintTime=TRUE,nb_t=nb,
                     tol=1e-3,maxIter=100,lowerBand=0.01,
                     upperBand=10)
    }
    obj
}
