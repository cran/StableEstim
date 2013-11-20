## Source file to design the main (S4) class Estim
## To be called at the top level of the package

## Validity check fct

validEstimObject <- function(object){
  ## check length of parameter estimation
  par <- object@par;
  if (length(par) ==4) ansp <- TRUE
  else ansp <- "Parameter of length different of 4"

  par0 <- object@par0;
  if (length(par0) ==4) ansp0 <- TRUE
  else ansp0 <- "Initial Parameter of length different of 4"

  ## check lenght covar
  vcov <- object@vcov;  
  if (ncol(vcov) ==4 && nrow(vcov)==4) anscov <- TRUE
  else anscov <- "covariance matrix of length different of 4x4"

  ## check length confint
  confint <- object@confint;  
  if (ncol(confint) ==2 && nrow(confint)==4) ansconfint <- TRUE
  else ansconfint <- "confidance intervall matrix of length different of 4x2"

  if (ansp==TRUE && ansp0==TRUE && anscov==TRUE && ansconfint==TRUE) res <- TRUE
  else if (is.character(ansp))   res <- ansp
  else if (is.character(ansp0))   res <- ansp0
  else if (is.character(anscov)) res <- anscov
  else if (is.character(ansconfint)) res <- ansconfint
  res
}

## Class Definition
setClass(Class="Estim",
         representation=representation(
             par        = "numeric",
             par0       = "numeric",
             vcov       = "matrix",
             confint    = "matrix", ## add an attr(confint,"level")<-level
             data       = "numeric",
             sampleSize = "numeric",
             others     = "list",
             duration   = "numeric",
             failure    = "numeric",
             method     = "character"),
         validity=validEstimObject
         )

## Init method

setMethod("initialize","Estim",
          function(.Object,par,par0,vcov,confint,
                   method,level,others,data,
                   duration,failure,...){
            ## handle missing
            if (missing(par))        par        <- numeric(4)
            if (missing(par0))       par0       <- numeric(4)
            if (missing(vcov))       vcov       <- matrix(nrow=4,ncol=4)
            if (missing(confint))    confint    <- matrix(nrow=4,ncol=2)
            if (missing(data))       data       <- numeric(100) 
            sampleSize <- length(data)
            if (missing(method))     method     <- "Default"
            if (missing(others))     others     <- list()
            if (missing(level))      level      <- 0
            if (missing(duration))   duration   <- 0
            if (missing(failure))    failure    <- 0
            
            ## set up names
            NameParamsObjects(par)
            NameParamsObjects(par0)
            NameParamsObjects(vcov)
            NameParamsObjects(confint);attr(confint,"level") <- level
            
            
            callNextMethod(.Object,par=par,par0=par0,vcov=vcov,
                           confint=confint,data=data,sampleSize=sampleSize,
                           method=method,others=others,duration=duration,
                           failure=failure,...)
        })

setMethod("show","Estim",
          function(object){
              cat("*** Class Estim, method Show *** \n")
              cat("** Method ** \n")
              print(object@method)
              cat("** Parameters Estimation ** \n")
              print(object@par)
              cat("** Covariance Matrix Estimation ** \n")
              print(object@vcov)
              cat("** Confidence interval Estimation ** \n")
              print(paste("Confidence level=",attributes(object@confint)$level))
              print(paste("data length=",object@sampleSize))
              print(object@confint)
              cat("** Estimation time ** \n")
              PrintDuration(object@duration)
              cat("** Estimation status ** \n")
              if(object@failure==0) cat("success")
              else cat("failure")
              cat("\n ******* End Show (Estim) ******* \n")
          }
          )
