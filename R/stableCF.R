ComplexCF <- function(t,theta,pm=0){
    alpha <- theta[1];beta <- theta[2]
    gamma <- theta[3];delta <- theta[4]
    CheckParametersRange(c(alpha,beta,gamma,delta))
    
    if ((pm !=0) & (pm !=1))
        stop(paste("pm= ",pm,"not taken into account. pm should be 0 or 1"))
    else ComputeComplexCF(t,alpha,beta,gamma,delta,pm)
}

jacobianComplexCF <- function(t,theta,pm=0){
    ComplexCFtoDeriv <- function(th)
        ComplexCF(t,th,pm)
    NumDeriv_jacobian(ComplexCFtoDeriv,theta)
}

ComputeComplexCF <- function(t,alpha,beta,gamma,delta,pm){
    arg <- t*delta
    scale <- -(gamma*abs(t))**alpha
    
    arg[t!=0] <-arg[t!=0]+ ComputeArgComplex(t,alpha,beta,gamma,delta,pm,scale)

    complex(real=exp(scale)*cos(arg),imaginary=exp(scale)*sin(arg))
}

ComputeArgComplex <- function(t,alpha,beta,gamma,delta,pm,scale){
    arg <- 0
    
    if (alpha == 1){
        if (gamma !=0){
            if (pm==0) arg <- scale[t!=0]*beta/pi2*sign(t[t!=0])*log(gamma*abs(t[t!=0]))
            else arg <- scale[t!=0]*beta/pi2*sign(t[t!=0])*log(abs(t[t!=0]))
        }
    }
    else {
        tanpa2 <- tan(pi2*alpha)
        if (gamma !=0){
            if (pm==0) arg <- scale[t!=0]*beta*tanpa2*sign(t[t!=0])*(abs(gamma*t[t!=0])**(1-alpha)-1)
            else arg <- -beta*tanpa2*sign(t[t!=0])*(scale[t!=0])
            }
    }
    arg   
}

## CHECK: need to be moved somewhere
CheckParametersRange <- function(theta){
    alpha <- theta[1]; beta <- theta[2];gamma <- theta[3]

    checkParams <- list(alpha=checkRange(alpha,0,2,"alpha"),
                     beta=checkRange(beta,-1,1,"beta"),
                     gamma=checkRange(gamma,0,ParamName="gamma")
                     )

    .printErr <- function(errList)
        if(!errList$bool) stop(errList$msg)

    lapply(checkParams,.printErr)
}

checkRange <- function(Parameter,min=-Inf,max=Inf,ParamName){
    if ((Parameter >= min) && (Parameter <= max)) return(list(bool=TRUE,msg="valid"))
    else return(list(bool=FALSE,msg=paste(ParamName, " = ",Parameter, " sould be in the interval ["
                                    ,min,max,"]")
                     )
                )  
}

## Empirical characteristic function evaluated at point t using data x
.ecf <- function(t,x)
{
    i <- complex(real=0,imaginary=1)
    fct <- function(tt) sum(exp(i*tt*x))/length(x)
    sapply(t,fct)
}
