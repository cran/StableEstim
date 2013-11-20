## We compute the C matrix as defined in the CGMM paper for length(x)=3
## we challenge the result against the integrate function (compute let by elt)

IntegrateRandomVectorsProduct.example <- function(){
    theta <- c(1.5,0.5,1,0)
    ## Integrands
    f_fct <- function(s,x){sapply(X=x,FUN=sampleComplexCFMoment,t=s,theta=theta)}
    f_bar_fct <- function(s,x){Conj(f_fct(s,x))}

    set.seed(345);X=rstable(3,1.5,0.5,1,0)
    ## Integration Params
    s_min=0;s_max=2
    numberIntegrationPoints=100
    randomIntegrationLaw="norm"
    
    Estim_Uniform <- IntegrateRandomVectorsProduct(f_fct,X,f_bar_fct,X,s_min,s_max,numberIntegrationPoints,
                                                  "Uniform",randomIntegrationLaw)
    Estim_Simpson <- IntegrateRandomVectorsProduct(f_fct,X,f_bar_fct,X,s_min,s_max,numberIntegrationPoints,
                                                  "Simpson",randomIntegrationLaw)
    ## Compute the result element by element using integrate
    mat <- matrix(0,3,3)
    Integrand_real <- function(s,i,j) Re(sampleComplexCFMoment(X[i],s,theta)*Conj(sampleComplexCFMoment(X[j],s,theta))*dnorm(s))
    Integrand_Im <- function(s,i,j) Im(sampleComplexCFMoment(X[i],s,theta)*Conj(sampleComplexCFMoment(X[j],s,theta))*dnorm(s))
    
    for (i in 1:3){
        for (j in 1:i){
            r <-integrate(f=Integrand_real,lower=s_min,upper=s_max,i=i,j=j)$value
            im <- integrate(f=Integrand_Im,lower=s_min,upper=s_max,i=i,j=j)$value
            mat[i,j] <- complex(real=r,imaginary=im)
            mat[j,i] <- Conj(mat[i,j])
        }
    }
    list(fct_unif=Estim_Uniform,
         fct_Simpson=Estim_Simpson,
         integrate=mat)
    ## User can check the accuracy of the Simpson scheme on this example
}
