RegularisedInverse.example <- function(){
    hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+")}

    K_h8 <- hilbert(8);
    r8 <- 1:8
    
    alphaReg_robust<- 1e-4
    Sa8_robust <- RegularisedSol(K_h8,alphaReg_robust,r8,"LF")
    
    alphaReg_accurate<- 1e-10
    Sa8_accurate <- RegularisedSol(K_h8,alphaReg_accurate,r8,"LF")
    
    ## when pre multiplied by K_h8 ,the expected solution is 1:8
    ## User can check the influence of the choice of alphaReg
    
    SolRobust <- K_h8 %*% Sa8_robust
    SolAccurate <- K_h8 %*% Sa8_accurate
    SolSolve <- K_h8 %*% solve(K_h8,r8)
    print("expected solution =1:8")
    list(robust=SolRobust,accurate=SolAccurate,solve=SolSolve)
}
