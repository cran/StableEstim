NumDeriv_jacobian <- function(fctToDeriv,WhereFctIsEvaluated,...){
    jacobian(fctToDeriv,WhereFctIsEvaluated,
             method="Richardson", method.args=list(),
             ...)
}
