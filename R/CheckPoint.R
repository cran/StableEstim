readCheckPoint <- function(AlphaBetaMatrix,Estimfct,nab,nSS,MCparam,...){
    method <- Estim_Des(Estimfct,...)
    fileName <- get_filename_checkPoint(AlphaBetaMatrix,nab,MCparam,method)

    if(!file.exists(fileName)){
        ## hearders line
        write(x="## ab;nab;sample;nSS;mc;MCparam",
              file=fileName,sep = "\n")
        ## initial vals
        ab <- 1;sample <- 1;mc <- 0
        write(x=paste("--",ab,nab,sample,nSS,mc,MCparam,sep=";"),
                  file=fileName,sep = "\n",append=TRUE)        
    }
    else {
        ## read current vals
        tab <- as.numeric(read.table(file=fileName,header=F,sep=";"))
        ab <- tab[2];sample <- tab[4];mc <- tab[6]
        n_ab <- tab[3];n_SS <-tab[5] ;mc_Param <- tab[7]
        stopifnot(nab==n_ab,nSS==n_SS,mc_Param==MCparam)
    }
    list(ab=ab,nab=nab,sample=sample,nSS=nSS,mc=mc,MCparam=MCparam)
}

writeCheckPoint <- function(AlphaBetaMatrix,Estimfct,ab,nab,
                            sample,nSS,mc,MCparam,...){
    method <- Estim_Des(Estimfct,...)
    fileName <- get_filename_checkPoint(AlphaBetaMatrix,nab,MCparam,method)
    
    line = readLines(fileName,-1)
    line[2]=paste("--",ab,nab,sample,nSS,mc,MCparam,sep=";")
    writeLines(line,fileName)
}

updateCheckPointValues <- function(CheckPointValues,MCparam,lS,nab){
    ab_start <- CheckPointValues$ab
    sample_start <- CheckPointValues$sample
    mc_start <- CheckPointValues$mc
    
    if (CheckPointValues$mc==MCparam){
        mc_start = 1
        if (CheckPointValues$sample==lS){
            sample_start = 1
            if (CheckPointValues$ab==nab)
                stop("Simulation finished already! check your output file")
            else ab_start = CheckPointValues$ab+1
        }
        else sample_start = CheckPointValues$sample +1
    }
    else mc_start = mc_start +1
    
    list(ab_start=ab_start,sample_start=sample_start,mc_start=mc_start)
}

deleteCheckPoint <- function(AlphaBetaMatrix,Estimfct,nab,
                             nSS,MCparam,...){
    method <- Estim_Des(Estimfct,...)
    fileName <- get_filename_checkPoint(AlphaBetaMatrix,nab,MCparam,method)

    unlink(x=fileName,force=TRUE)
}

get_filename_checkPoint <- function(alphaBetaMat,nab,MCparam,method){
    
    MC <- paste(paste("alpha0=",alphaBetaMat[1,1],sep=""),
                paste("beta0=",alphaBetaMat[1,2],sep=""),
                paste("alphan=",alphaBetaMat[nab,1],sep=""),
                paste("betan=",alphaBetaMat[nab,2],sep=""),
                paste("MCparam",MCparam,sep=""),
                sep="_"
                  )
    fileName <- paste(MC,method,"_CHECKPOINT.txt",sep="")
    fileName
}
