test.ComputeWeightingMatrix <- function(){
    t <- seq(0.1,1.5,length.out=5)
    theta <- c(1.4,0.4,1,0)
    set.seed(345); x <- rstable(2500,1.,0.4)

    test.OptAsym.case(t,theta,x)
    test.DataVar.case(t,theta,x)
    test.id.case(t,theta,x)
}

test.OptAsym.case <- function(t,theta,x){
    Kestim <- ComputeWeightingMatrix(t,theta,x,"OptAsym")
    trKestim <- sum(diag(Kestim))
    trRef <- 3.207883

    expect_almost_equal(trKestim,trRef)
}

test.DataVar.case <- function(t,theta,x){
    Kestim <- ComputeWeightingMatrix(t,theta,x,"DataVar")
    trKestim <- sum(diag(Kestim))
    trRef <- 2.79743

    expect_almost_equal(trKestim,trRef)
}

test.id.case <- function(t,theta,x){
    Kestim <- ComputeWeightingMatrix(t,theta,x,"Id")
    trKestim <- sum(diag(Kestim))
    expect_equal(trKestim,2*length(t))
}

ComputeWeightingMatrix <- function(t,theta,x,WeightingMatrix,pm=0,...){
    switch(WeightingMatrix,
           "OptAsym"={K <- asymVarRealCFMoment(t=t,theta=theta,pm=pm)},
           "DataVar"={K <- DataVarRealCFMoment(t=t,theta=theta,x=x,pm=pm)},
           "Id"={K <- diag(nrow=2*length(t),ncol=2*length(t))},
           "provided"={
               if(!is.null(Kt <- list(...)$ProvidedWeightingMatrix)) K <- Kt
               else  stop("You need to provide a Weighting matrix")}
           )
    K
}

## Km=E[g(x,theta)g(x,theta)']=var[g(V,theta)] g is CFbased
## formula due to Qin and Lawless (1994) (their Lemma 1 and Theorem 1) 
## and recalled by Owada (2006)

test.asymVarRealCFMoment <- function(){
    t <- seq(0.1,1.5,length.out=5)
    theta <- c(1.4,0.4,1,0)
    
    Kestim <- asymVarRealCFMoment(t,theta)
    trKestim <- sum(diag(Kestim))
    trRef <- 3.207883

    expect_almost_equal(trKestim,trRef)
}

asymVarRealCFMoment<- function(t,theta,pm=0){
    m    <- length(t)
    res  <- matrix(0,ncol=2*m,nrow=2*m)
    
    tiPlustj <- crossSum(t,t) 
    tiMenostj <- crossSum(t,-t)
    phi  <- ComplexCF(t,theta,pm) 
    phi_tiPlus_tj <- sapply(tiPlustj,ComplexCF,theta,pm)
    phi_tiMenos_tj <- sapply(tiMenostj,ComplexCF,theta,pm)
    
    res[1:m        ,1:m        ] <- (0.5*(Re(phi_tiPlus_tj)+ Re( phi_tiMenos_tj))
                                     - Re(phi) %*% t(Re(phi))) # block(1m,1m)
    res[1:m        ,(m+1):(2*m)] <- (0.5*(Im(phi_tiPlus_tj)- Im( phi_tiMenos_tj))
                                     - Re(phi) %*% t(Im(phi))) # block(1m,m+1:2m)
    res[(m+1):(2*m),1:m        ] <- (0.5*(Im(phi_tiPlus_tj)+ Im( phi_tiMenos_tj))
                                     - Im(phi) %*% t(Re(phi))) # block(m+1:2m,1m)
    res[(m+1):(2*m),(m+1):(2*m)] <-(-0.5*(Re(phi_tiPlus_tj)- Re( phi_tiMenos_tj))
                                    - Im(phi) %*% t(Im(phi))) # block(m+1:2m,m+1:2m)
    res
}

DataVarRealCFMoment<- function(t,theta,x,pm=0){
    gt <- DataMatrixRealCFMomentCondition(t=t,theta=theta,x=x,pm=pm)
    var(gt)
}

test.DataVarRealCFMoment <- function(){
    t <- seq(0.1,1.5,length.out=5)
    theta <- c(1.4,0.4,1,0)
    set.seed(345); x <- rstable(2500,1.,0.4)
    
    Kestim <- DataVarRealCFMoment(t,theta,x)
    trKestim <- sum(diag(Kestim))
    trRef <- 2.79743
    
    expect_almost_equal(trKestim,trRef)

}
