sampleComplexCFMoment <- function(x,t,theta,pm=0){
    ecf <- function(tt) mean(exp(complex(imaginary=tt*x)))
    phiXn <- sapply(t,ecf)
    phiTheta <- ComplexCF(t,theta,pm)
    return(phiXn-phiTheta)
}

## res=c(Real,Im) hence size=2*size(t)
sampleRealCFMoment <- function(x,t,theta,pm=0){
    ComplexRes <- sampleComplexCFMoment(x,t,theta,pm)
    return(c(Re(ComplexRes),Im(ComplexRes)))
}

jacobianSampleComplexCFMoment<- function(t,theta,pm=0){
    -jacobianComplexCF(t,theta,pm)
}

jacobianSampleRealCFMoment<- function(t,theta,pm=0){
    jac <- jacobianSampleComplexCFMoment(t,theta,pm)
    apply(jac,2,function(x) c(Re(x),Im(x)))
}

DataMatrixRealCFMomentCondition <- function(t,theta,x,pm=0){
    x <- matrix(c(x),ncol=1)
    x_comp <- x%*%matrix(t,nrow=1)
    x_comp <- matrix(complex(imaginary=x_comp),ncol=length(t))
    emp_car <- exp(x_comp)
    the_car <- ComplexCF(t,theta,pm)
    gt <- t(t(emp_car) - the_car)
    gt <- cbind(Re(gt),Im(gt))
    return(gt)
}
