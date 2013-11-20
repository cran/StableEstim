## We compute the regularized solution of the original problem K \phi =r
## the (non-regularized solution) phi=\sum_{j=1}^\infty \frac{<r,\\ksi_j>}{\lambda_j} \phi_j
## where \lambda_j is the singular values
##       \ksi_j left singular vector
##       \phi_j right singular vector
## The regularized solution is given by R_{\alpha,r}= \sum_{j=1}^\infty q(\alpha,\lambda_j)\frac{<r,\\ksi_j>}{\lambda_j} \phi_j
## where q(\alpha,\lambda_j) is one of the regularization scheme described in the paper
## LF is the Landweber-Fridmann and cut-off the spectral cut-off

RegularisedSol <- function(Kn,alphaReg,r,regularization=c("Tikhonov","LF","cut-off"),...){
    regularization <- match.arg(regularization)
    return(ComputeRegularizedSol(Kn=Kn,alpha=alphaReg,r=r,regularization=regularization,...))

}

ComputeRegularizedSol <- function(Kn,alpha,r,regularization=c("Tikhonov","LF","cut-off"),...){
    regularization<- match.arg(regularization)

    singularValuesSystem <- getSingularValueDecomposition(Kn)
    lambda <- singularValuesSystem$lambda
    phi <- singularValuesSystem$phi
    ksi <- singularValuesSystem$ksi

    qAlphaLambda <- ComputeqAlphaLambda(lambda=lambda,alpha=alpha,
                                        regularization=regularization,
                                        ...,Kn=Kn)
    SP_rbyKsi <- crossprod(r,ksi)
    qAlphaSP <- t(qAlphaLambda/lambda) *SP_rbyKsi
    return(rowSums(matrix(data=qAlphaSP,ncol=ncol(phi),nrow=nrow(phi),byrow=TRUE)*phi))
}

## if Kn is a symmetric matrix, we use eigenvalues analysis
getSingularValueDecomposition <- function(Kn){
    if (isSymmetric(Kn)){
        SingularValuesDecomposition <- eigen(x=Kn, symmetric=TRUE)
        phi <- SingularValuesDecomposition$vectors
        ksi <- SingularValuesDecomposition$vectors
        lambda <- SingularValuesDecomposition$values
    }
    else {
        SingularValuesDecomposition <- svd(Kn)
        phi <- SingularValuesDecomposition$v
        ksi <- SingularValuesDecomposition$u
        lambda <- SingularValuesDecomposition$d
    }
    return(list(lambda=lambda,phi=phi,ksi=ksi))
}

ComputeqAlphaLambda <- function(lambda,alpha,regularization,...){
    switch(regularization,
           "Tikhonov"={
               lambda2 <- lambda**2
               return(lambda2/(lambda2+alpha))
           },
           "LF"={
             c <- checkValidConstantC(...)
             invAlpha <- floor((1/alpha))
             return(1-(1-c*(lambda**2))**invAlpha)
           },
           "cut-off"={
               sqrtAlpha <- sqrt(alpha)
               TestCutOff <- function(x,a){ifelse(x<sqrt(a),x**2/a, 1)}
               return(sapply(lambda,TestCutOff,a=alpha))
           },
           stop(paste(regularization," method not taken into account"))
           )
}


checkValidConstantC <- function(...){
    eps <- 1e-6
    args <- list(...);
    c <- args$c
    K <- args$Kn; if (is.null(K)) stop("you need to provide Kn for LF Regularization")
    InvNormK2 <- 1/(norm(K)**2)
    if ((is.null(c)) || (c < 0) || (c > InvNormK2) ){
        c <- max(eps,InvNormK2-eps)}
    else
        c <- c
    return(c)
    }

ComputeCutOffInverse <- function(X,alphaReg=0.001){
    s <- getSingularValueDecomposition(X)
    index <-(abs(s$lambda) < alphaReg)
    if (any(index)){
        lambda <- s$lambda; lambda[index] <- alphaReg
        D <- diag(lambda)
        Invmat <- s$ksi %*% D %*% t(s$phi)
    }
    else{
        qx <- qr(X)
        Invmat <-solve.qr(qx)
    }
    Invmat
}
