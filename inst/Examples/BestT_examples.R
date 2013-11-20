#=================================================================================
#------------------------------- Best T ------------------------------------------
#=================================================================================
ComputeBest_t.example <- function(){
    AlphaBetaMatrix <- matrix(c(0.7,1.5,0,0.5),ncol=2)
    nbts <- seq(20,50,10)
    
    ComputeBest_t(AlphaBetaMatrix=AlphaBetaMatrix,
                  nb_ts=nbts)
}

#=================================================================================
#------------------------------- Best tau ------------------------------------------
#=================================================================================
ComputeBest_tau.example <- function(){
    AlphaBetaMatrix <- matrix(c(0.7,1.5,0,0.5),ncol=2)
    nbts <- seq(20,50,10)
    
    ComputeBest_tau(AlphaBetaMatrix=AlphaBetaMatrix,
                    nb_ts=nbts,tScheme="uniformOpt",
                    Constrained=TRUE)
    

}
