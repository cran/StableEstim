##==========================================================================================
##--------------------------------- MC Simulation ------------------------------------------
##==========================================================================================
MCSim.example <- function(){ #CHECK: need to export some function
    ## Method Params
    algo="2SGMM"
    alphaReg=0.05
    regularization="Tikhonov"
    WeightingMatrix="DataVar"
    t_scheme="uniformOpt"
    pm=0

    ## Simulation Params
    AlphaBetaMatrix <- matrix(c(0.7,1.5,0,0.5),ncol=2)
    SampleSizes <- c(300,600,1000)
    MCparam <- 10
    Estimfct <- "GMM"

    Estim_Simulation(AlphaBetaMatrix=AlphaBetaMatrix,MCparam=MCparam,SampleSizes=SampleSizes,
                     Estimfct=Estimfct,HandleError=FALSE,saveOutput=TRUE,StatSummary=TRUE,
                     algo=algo,alphaReg=alphaReg,
                     regularization=regularization,WeightingMatrix=WeightingMatrix,
                     t_scheme=t_scheme,pm=pm,nb_t=42,Constrained=FALSE,
                     tol=1e-3,maxIter=100,lowerBand=0.01,
                     upperBand=10)

    ## in order to check the checkPointing process, you can kill the procedure during
    ## running time and run it again and check how the process reacts
}



##==========================================================================================
##--------------------------------- File Manipulation --------------------------------------
##==========================================================================================
fileManip.example <- function(){
    file1="file1.csv" # No header
    file2="file2.csv" # header
    files <- c(file2,file1)
    
    ## ----------------Concatenates files in a file called fileSum.csv -------
    ConcatFiles(files=files,outfile="fileSum.csv",headers_=c(TRUE,FALSE),DeleteIfExists=TRUE)
    
    ## ---------------Create stat objects for concat file ------------------------------------
    fileObj <- ComputeStatObjectFromFiles("fileSum.csv")
    fileObj
}

##==========================================================================================
##--------------------------------- Tex Summary --------------------------------------
##==========================================================================================
texSummaryFromObj.example <- function(){
    obj <- fileManip.example()
    texObj <- TexSummary(obj,par_index=1:4,digits=c(3,3,3,4,4,4))

    list(Obj=obj,tex=texObj)
}

texSummaryFromFile.example <- function(){
    file1="test1.csv" # No header
    file2="test2.csv" # header
    files <- c(file2,file1)

    
    TexSummary(files=files,par_index=1:4,digits=c(3,3,3,4,4,4),
               headers_=TRUE)
}


