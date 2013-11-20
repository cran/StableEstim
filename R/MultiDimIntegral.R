## We want to compute vectorial Integrals of the form C=\int_a^b f_s(X) g_s(Y) \pi(s) ds
## where f_s \in \mathbb{R}^n (X \in \mathbb{R}^m) with n \gep 1 (m \gep 1) (same works for g and Y)
## the elements of C are given by c_ij=\int_a^b f_s(X_i) g_s(Y_j) \pi(s) ds
## We implement 2 methods using different integration scheme:
## 1) Trapezian scheme ; 2) Simpson scheme
## The previous procedure are sorted by a speed decresing criteria(vs incresing accuracy)
## We expect f ( and g) to be functions that return a matrix of the size length(s)xlength(X)
## i.e f=(f_si(x_j))_{1 \leq i \leq length(s) , 1 \leq j \leq length(X)}  

IntegrateRandomVectorsProduct <- function(f_fct,X,g_fct,Y,s_min,s_max,
                                          subdivisions=50,
                                          method=c("Uniform","Simpson"),
                                          randomIntegrationLaw=c("norm","unif"),...){
    randomIntegrationLaw <- match.arg(randomIntegrationLaw)
    method <- match.arg(method)

    intBlock<- function(sequence){getWeightedVectorProductForSpecifiedSequnce(f_fct=f_fct,X=X,g_fct=g_fct,Y=Y,
                                                                              s_min=s_min,s_max=s_max,
                                                                              subdivisions=subdivisions,
                                                                              randomIntegrationLaw=randomIntegrationLaw,
                                                                              sequence=sequence,...)}
    switch(method,
           "Uniform"={
               int <-intBlock("uniform")
           },
           "Simpson"={
               int2j <-intBlock("2j")
               int2j1<-intBlock("2j-1")
               intmin<- intBlock("min")
               intmax <- intBlock("max")
               int <- (intmin+ 2*int2j + 4*int2j1 + intmax)/3
           },
           stop(paste("Integration method",method," not taken into account"))
           )
    int
}

getWeightedVectorProductForSpecifiedSequnce <- function(f_fct,X,g_fct,Y,s_min,s_max,
                                                        subdivisions,
                                                        randomIntegrationLaw,
                                                        sequence=c("uniform","2j","2j-1","min","max"),...){
    sequence <- match.arg(sequence)
    ns <- subdivisions
    h <- (s_max-s_min)/ns
    
    s <- switch(sequence,
                "uniform"={seq(from=s_min,to=s_max,length.out=(ns+1))}, ## nbr int points= nbr subintervalss +1
                "2j"={seq(s_min+2*h,s_min+(ns-2)*h,2*h)},
                "2j-1"={seq(s_min+h,s_min+(ns-1)*h,2*h)},
                "min"={s_min},
                "max"={s_max},
                stop(paste("sequence",sequence," not taken into account"))
                )
    
    f <- matrix(data=f_fct(s,X),nrow=length(s)) # matrix ns x n_X or ns x n_X
    g <- matrix(data=g_fct(s,Y),nrow=length(s)) # matrix ns x n_Y or ns x n_Y

    
    randomWeight_f <- sqrt(as.vector(getRandomWeight(IntegrationPoints=s,lengthOutput=nrow(f),
                                                     randomIntegrationLaw=randomIntegrationLaw,
                                                     s_min=s_min,s_max=s_max,...)*h
                                     )
                           )
    randomWeight_g <- sqrt(as.vector(getRandomWeight(IntegrationPoints=s,lengthOutput=nrow(g),
                                                     randomIntegrationLaw=randomIntegrationLaw,
                                                     s_min=s_min,s_max=s_max,...)*h
                                     )
                           )
                                     

    fByRandomWeight <- f*randomWeight_f
    gByRandomWeight <- g*randomWeight_g
    
    return(crossprod(fByRandomWeight,gByRandomWeight))
}

getRandomWeight <- function(IntegrationPoints,lengthOutput,
                            randomIntegrationLaw,s_min,s_max,...){
    args <- list(...)
    switch(randomIntegrationLaw,
           "unif"={
               lower <- s_min
               upper <- s_max
               res <- dunif(x=IntegrationPoints,min=lower,max=upper)  
           },
           "norm"={
               mu <- ifelse(is.null(args$mu),0,args$mu)
               sigma <- ifelse(is.null(args$sigma),1,args$sigma)
               res <- dnorm(x=IntegrationPoints,mean=mu,sd=sigma)
           },
           stop(paste("random Integration Law",randomIntegrationLaw," not taken into account"))
           )
    
    matrix(data=res,nrow=lengthOutput,ncol=1)
}


