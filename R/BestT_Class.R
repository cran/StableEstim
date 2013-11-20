validEstimObject_t<- function(object){
  ## check length of parameter estimation
  par <- object@theta;
  if (length(par) ==4) ansp <- TRUE
  else ansp <- "Parameter of length different of 4"

  ansp
}

## Class Definition
setClass("Best_t",
         representation(theta   = "vector",
                        nbt     = "vector",
                        tvec    = "list",
                        detVal  = "vector",
                        convcode= "numeric"),
         validity=validEstimObject_t)

## Init method
setMethod("initialize","Best_t",
          function(.Object,theta,nbt,
                   tvec,detVal,convcode,...){
                   
            ## handle missing
            if (missing(theta))     theta     <- vector(mode="numeric",length=4)
            if (missing(nbt))       nbt       <- seq(10,100,10)
            if (missing(tvec))      tvec      <- lapply(nbt,function(n) numeric(n))
            if (missing(detVal))    detVal    <- vector(mode="numeric",length=length(nbt))
            if (missing(convcode))  convcode  <- 0
            
            callNextMethod(.Object,theta=theta,nbt=nbt,
                           tvec=tvec,detVal=detVal,
                           convcode=convcode,...)
        })

setMethod("+", signature(e1 = "Best_t", e2 = "Best_t"), function (e1, e2){
    if (all(e1@theta==e2@theta)){
        res <- new("Best_t",theta=e1@theta,
                   nbt=sort(c(e1@nbt,e2@nbt)))
        
        orderIndex <- order(c(e1@nbt,e2@nbt))
        
        for (i in 1:length(res@nbt)){
            Global_ind <- orderIndex[i]
            if (Global_ind > length(e1@nbt)){
                ind <- Global_ind - length(e1@nbt)
                res@tvec[[i]] <- e2@tvec[[ind]]
                res@detVal[i] <- e2@detVal[ind]
            }
            else{
                ind <- Global_ind
                res@tvec[[i]] <- e1@tvec[[ind]]
                res@detVal[i] <- e1@detVal[ind]
            }        
        }
        
        names(res@detVal) <- res@nbt
        res@convcode <- max(e1@convcode,e2@convcode)
        
        return(res)
    }
    else stop("cannot sum to vectors corresponding to different values of theta")
}
          )

setMethod("show","Best_t",
          function(object){
              cat("*** Class Best_t, method Show *** \n")
              cat("** theta ** \n")
              print(object@theta)
              cat("** number of points ** \n")
              print(object@nbt)
              cat("** determinant value  ** \n")
              print(object@detVal)
              cat("** Convergence code ** \n")
              print(object@convcode)
              cat("\n ******* End Show (Estim) ******* \n")
          }
          )
