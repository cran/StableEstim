## each files contains a unique set of Parametes: one and only one file per set of parameters
## if you have more than one file for the same set of True value, use ConcatFiles first
## We assume the same set of sample sizes to be used in all the files:
## those sample sizes are read from file n=readSizeFrom
ComputeStatObjectFromFiles <- function(files,sep_=",",FctsToApply=StatFcts,
                                       headers_=TRUE,readSizeFrom=1,
                                       CheckMat=TRUE,tolFailCheck=tolFailure,
                                       MCparam=1000,...){
    nab <- length(files)

    headers <- headers_
    sep <- sep_
    if (nab > 1){
        if (length(sep_)!=nab) sep <- rep(sep_,nab)
        if (length(headers_)!=nab) headers <- rep(headers_,nab)
    }
    
    SampleSizes <- getSampleSizesFromFiles(files[readSizeFrom],sep)
    
    lS <- length(SampleSizes)
    nRowOutput <- nab*lS
    indexStatOutput <- 1
    StatOutputLength <- length(FctsToApply)+ 5 # 5=Sample size+alphaT+betaT+Total failure+average time

    StatOutput <- list(alpha=matrix(data=NA,ncol=StatOutputLength,nrow=nRowOutput),
                       beta =matrix(data=NA,ncol=StatOutputLength,nrow=nRowOutput),
                       gamma=matrix(data=NA,ncol=StatOutputLength,nrow=nRowOutput),
                       delta=matrix(data=NA,ncol=StatOutputLength,nrow=nRowOutput))
    
    for (ab in 1:nab){
        data <- as.matrix(read.csv(file=files[ab],sep=sep[ab],header=headers[ab]))
        
        res <- ComputeStatOutput(EstimOutput=data,
                                 FctsToApply=FctsToApply,
                                 SampleSizes=SampleSizes,
                                 CheckMat=CheckMat,
                                 tolFailCheck=tolFailCheck,
                                 MCparam=MCparam,...)
        IndexSec <- seq(indexStatOutput,indexStatOutput+(lS-1),1)
        StatOutput$alpha[IndexSec,] <- res$alpha
        StatOutput$beta [IndexSec,] <- res$beta
        StatOutput$gamma[IndexSec,] <- res$gamma
        StatOutput$delta[IndexSec,] <- res$delta
        indexStatOutput <- indexStatOutput +lS
    }
    NameStatOutput(FctsToApply,StatOutput)
}

## Concatenate files output in one single file
## checks if the set of true values is the same
ConcatFiles <- function(files,sep_=",",outfile,headers_=TRUE,
                        DeleteIfExists=TRUE){
    nab <- length(files)
    header_flag=FALSE
    Exists=file.exists(outfile)
    
    if (DeleteIfExists && Exists) unlink(x=outfile,force=TRUE)
    
    headers <- headers_
    sep <- sep_
    if (nab > 1){
        if (length(sep_)!=nab) sep <- rep(sep_,nab)
        if (length(headers_)!=nab) headers <- rep(headers_,nab)
    }

    stopifnot(CheckTrueParInFiles(files,sep,headers)$answer)
    file.create(outfile)
    writeLines("",outfile)
    
    for (ab in 1:nab){
        data <- as.matrix(read.csv(file=files[ab],sep=sep[ab],header=headers[ab]))
        if(headers[ab] && !header_flag){
            header_flag=TRUE
            ## write headers in first line
            names <- colnames(data)
            if (is.null(names)) stop("There is no headers! Check headers_ arg!!")
            line = readLines(outfile,-1)
            l1=names[1]
            for (i in 2:length(names)) l1 <- paste(l1,names[i],sep=sep[ab])
            line[1]=l1
            writeLines(line,outfile)
            ## append file
            write.table(x=data,file=outfile,sep=sep[ab],
                        col.names=FALSE,row.names=FALSE,
                        append=TRUE)
        }
        else {
            write.table(x=data,file=outfile,sep=sep[ab],
                        col.names=FALSE,row.names=FALSE,
                        append=TRUE)
        }
    }   
}

CheckTrueParInFiles <- function(files,sep,headers){
    nab <- length(files)
    TrueParMat <- matrix(nrow=nab,ncol=2)
    
    for (ab in 1:nab){
        data <- as.matrix(read.csv(file=files[ab],
                                   sep=sep[ab],
                                   header=headers[ab]
                                   )
                          )
        TrueParCheck <- checkTrueParInMatrix(data,files[ab])
        TrueParMat[ab,] <- c(TrueParCheck$alpha,
                             TrueParCheck$beta)
    }
    checkTrueParInMatrix(TrueParMat)
}

checkTrueParInMatrix <- function(data,file=NULL){
    alpha=data[1,1]
    beta=data[1,2]

    answers <- apply(data,1,
                     function(x){ifelse((x[1]==alpha) && (x[2]==beta),TRUE,FALSE)}
                     )
    answer <- ifelse(any(answers==FALSE),FALSE,TRUE)
    if (!is.null(file) && !answer)
        stop(paste("file=",file,"contains different values of alpha and beta",sep=" ")) 
    list(alpha=alpha,beta=beta,answer=answer)
}

## produce a list of objects from class "Latex" using either:
## 1) the object produced Estim_Simulation
## 2) the file saved by Estim_Simulation
## index could be a number in 1:4, a sequence within 1:4 or
## one or more value in c("alpha","beta","gamma","delta")

TexSummary <- function(obj,files=NULL,sep_=",",FctsToApply=StatFcts,
                       caption="Statistical Summary",
                       label='Simtab',digits=3,par_index=1,MCparam=1000,...){
    output=list()
    if (missing(obj)){
        if(is.null(files)) stop("You need to input a valid Object or specify a file-name to parse !!")
        else obj <- ComputeStatObjectFromFiles(files=files,sep_=sep_,
                                               FctsToApply=FctsToApply,
                                               MCparam=MCparam,...)
    }
    format <- get_format(obj[[par_index[1]]],digits,FctsToApply)
    
    for (i in 1:length(par_index)){
        mat <- obj[[par_index[i]]]
        tobj <- xtable(x=mat,
                       caption=caption[i],
                       label=label[i])
        
        digits(tobj) <- format$digits 
        align(tobj)  <- format$align
        display(tobj) <- format$display
        output[[par_index[i]]] <- toLatex(tobj,include.rownames=FALSE)
    }
    output
}

get_format <- function(mat,digits,FctsToApply){
    nf <- length(FctsToApply)
    nmat <- ncol(mat)
    stopifnot((nmat-nf)==5)

    rdigits <- get_digits(digits,nf,mat)
    align <- get_align(nmat,nf)
    display <- get_display(nmat,nf)
    
    list(digits=rdigits,align=align,display=display)
}

get_align <- function(nmat,nf){
    res <-  character(nmat+6) #4=1st col  + 3 |
    res[1:7] <- c("c","|","c","c","|","c","|")
    res[8:(8+nf-1)] <- rep("c",nf)
    res[(8+nf):(nmat+6)] <- c("|","c","c","|")
    
    as.character(res)
}

get_display <- function(nmat,nf){
    res <-  character(nmat+1)
    res[c(4,nmat)] <- "d"
    res[c(1:3,5:(5+nf-1),nmat+1)] <- "f"
    
    as.character(res)
}

get_digits <- function(digits,nf,mat){
    nmat <- ncol(mat)
        
    if (is.matrix(digits)){
        if (length(mat) != length(digits)) stop("check digits matrix size !!")
        else {res <- matrix(ncol=ncol(digits)+1,nrow=nrow(digits))
              res[,1]=0
          }
    }
    else {
        if (length(digits)==nf){res <- c(1,3,3,4,digits,3,4)}
        else if (length(digits)==nmat){res <- c(1,digits)}
        else if (length(digits)==1){res <- c(1,2,2,0,rep(digits,nf),0,2)}
        else stop("digits size is not correct !!")
    }
    res
}

saveFile <- function(Output,alphaT,betaT,MCparam){
    method <- Output$file
    fileName <- get_filename(alphaT,betaT,MCparam,method)
    
    v <- write.table(x=Output$outputMat,file=fileName,
                     sep = ",",row.names = FALSE,
                     col.names=TRUE)
}

initOutputFile <- function(alphaT,betaT,MCparam,Estimfct,...){
    method <- Estim_Des(Estimfct,...)
    fileName <- get_filename(alphaT,betaT,MCparam,method)

    if(!file.exists(fileName)){
        write(x=paste("alphaT","betaT","data size","seed","alphaE",
                  "betaE","gammaE","deltaE","failure","time",sep=","),
              file=fileName,sep = "\n")
    }
}

updateOutputFile <- function(alphaT,betaT,MCparam,Output){
    method <- Output$file
    fileName <- get_filename(alphaT,betaT,MCparam,method)
    
    write(x=paste(as.character(Output$outputMat),collapse=","),
          file=fileName,sep="\n",append=TRUE)
}

get_filename <- function(alphaT,betaT,MCparam,method,extension=".csv"){
    MC <- paste(paste("alphaT=",alphaT,sep=""),
                paste("betaT=",betaT,sep=""),
                paste("MCparam",MCparam,sep=""),
                sep="_"
                  )
    fileName <- paste(MC,method,extension,sep="")
    fileName
}

getSampleSizesFromFiles <- function(file,sepa=","){
    data <- as.matrix(read.csv(file=file,sep=sepa))
    if (!is.numeric(data))
        stop(paste("Check the file",file, "! Output is not numeric when sep=",sepa, " is used."),sep="")
    unique(data[,3])
}

NameStatOutput <- function(FctsToApply,StatOutput) {
    Names <- c("alpha","beta","n",names(FctsToApply),"failure","time")
    lapply(StatOutput,function(x) {colnames(x) <- Names;return(x)})
}
