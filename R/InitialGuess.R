IGParametersEstim<-function(x,pm=0,...) {
    par <- KogonParametersEstim(x,pm)
    alpha <- par[1];beta <- par[2];gam <- par[3]
    if (alpha==2) alpha <- alpha-eps
    if (alpha==0) alpha <- eps
    if (beta==-1) beta <- -1+eps
    if (beta==1) beta <- 1-eps
    if (gam==0) gam <- eps
    c(alpha,beta,gam,par[4])
}

McCullochParametersEstim<-function(x){
        tr <- tryCatch(.qStableFit(x, doplot=FALSE),
                       error=function(e)e) 
        err <- inherits(tr, "error")
        if (!err) {
            ans <- tr
            alpha <- ifelse(is.na(ans@fit$estimate[1]),0.5,ans@fit$estimate[1])
            beta  <- ifelse(is.na(ans@fit$estimate[2]),0,ans@fit$estimate[2])
            gamma <- ifelse(is.na(ans@fit$estimate[3]),1,ans@fit$estimate[3])
            delta <- ifelse(is.na(ans@fit$estimate[4]),0,ans@fit$estimate[4])
        }
        else{
            alpha <- 0.5
            beta <- 0
            gamma <- 1;delta <- 0
        }
        c(alpha,beta,gamma,delta)
    }

KogonParametersEstim<-function(x,pm=0){
    t <- seq(from=0.1,to=1,length.out=10)
    Par <- McCullochParametersEstim(x)
    alpha <- Par[1];beta <- Par[2]
    gamma <- Par[3];delta <- Par[4]

    ## Shift data to fit the regression requirement
    if (pm==0) x <- x+beta*gamma*tan(pi2*alpha)
    
    s <- (x-delta)/gamma
    PhiHat <- .ecf(t,s)
    alpha <- EstimateAlpha(PhiHat,t)
    beta <- EstimateBeta(s,t,alpha)

    c(alpha=alpha,beta=beta,gamma,delta)
}

EstimateAlpha <- function(PhiHat,t){
    Y1 <- log(-log(Mod(PhiHat)))
    X1 <- log(Mod(t))
    reg1 <- lm(Y1~X1)
    alpha <- as.numeric((reg1[[1]])[2])
    alpha=max(alpha,0)
    alpha=min(alpha,2)
    alpha
}

EstimateBeta <- function(s,t,alpha){
    Y2 <- .get.arctan(s,t)
    X2 <- (abs(t)**alpha)*(t/abs(t))*tan(pi2*alpha)
    reg2 <- lm(Y2~X2)
    beta <- as.numeric((reg2[[1]])[2])
    beta=min(beta,1);beta=max(beta,-1)
    beta
}
