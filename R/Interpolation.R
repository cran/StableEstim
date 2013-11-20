## xTab should be monotically increasing or decreasing
## Given xTab,m <=n and xVal, find an integer l0 such that
## xVal is centered among x_lo,...,x(l0+M-1)
## m = l0 + floor((m-2)/2)
## l0 cannot be less than 0
## l0+M-1 cannot be greater thab N-1.

linearInterp_1D <- function(xTab,yTab,xVal,decreasing=FALSE){
    ## Init 
    xTab <- as.numeric(xTab)
    n <- length(xTab)
    or <- order(xTab,decreasing=decreasing)
    SxTab <- xTab[or]
    SyTab <- yTab[or]
    isav <- 1
    dj <- min(1,n**0.25)
    cor <- 0

    ## Linear Interp specific
    m <- 2
    l0 <- ifelse(cor,
                 MixedLocalize(SxTab,xVal,m,dj,isav,n),
                 BisectionLocalize(SxTab,xVal,m,dj,isav,n)
                 )

    if (SxTab[l0]==SxTab[l0+1]) return (SyTab[l0])
    else return(SyTab[l0]+ ((xVal-SxTab[l0])/(SxTab[l0+1]-SxTab[l0]))*(SyTab[l0+1]-SyTab[l0]))
}

BilinInterp_2D <- function(xTab,yTab,zMat,xVal,yVal,decreasing=FALSE){
    if (xVal > max(xTab) || yVal > max(yTab)) return(NA)
    if (xVal < min(xTab) || yVal < min(yTab)) return(NA)
    
    N <- length(as.numeric(xTab))
    M <- length(as.numeric(yTab))
    stopifnot(length(zMat)==N*M)

    zMat <- matrix(zMat,nrow=N,ncol=M)
    ## get i
    xTab <- as.numeric(xTab)
    isav <- 1
    dj <- min(1,N**0.25)
    i <- BisectionLocalize(xTab=xTab,xVal=xVal,m=2,dj=dj,isav=isav,n=N)

    ## get j
    yTab <- as.numeric(yTab)
    isav <- 1
    dj <- min(1,M**0.25)
    j <- BisectionLocalize(xTab=yTab,xVal=yVal,m=2,dj=dj,isav=isav,n=M)

    ## interpolate
    t <- (xVal-xTab[i])/(xTab[i+1]-xTab[i])
    u <- (yVal-yTab[j])/(yTab[j+1]-yTab[j])

    return((1-t)*(1-u)*zMat[i,j]+t*(1-u)*zMat[i+1,j]+(1-t)*u*zMat[i,j+1]+t*u*zMat[i+1,j+1])
}

BisectionLocalize <- function(xTab,xVal,m,dj,isav=1,n=length(xTab))
    {
        stopifnot(n >= 2, m >= 2, m <=n)

        ascnd <- (xTab[n] >= xTab[1])
        il <- 1;
        iu <- n
        while(iu-il > 1){
            im=(iu+il) %/% 2
            if ((xVal >= xTab[im])== ascnd)
                il=im
            else
                iu=im
        }
        cor <- ifelse(abs(il - isav) > dj, 0,1)
        isav=il
        return(max(0,min(n-m,il-((m-2) %/% 2))))
    }

MixedLocalize <- function(xTab,xVal,m,dj,isav=1,n=length(xTab))
    {
        il <- isav
        inc <- 1

        stopifnot(n >= 2, m >= 2, m <=n)
        ascnd <- (xTab[n] >= xTab[1])

        if (il < 1 || il > n){
            il=1
            iu=n
        }
        else {
            if ((xVal >=xTab[il]) == ascnd) {
                for (index in 1:n){
                    iu=il+inc
                    if (iu >= n) {iu=n; break}
                    else if ((xVal < xTab[iu]) == ascnd) break
                    else {
                        il=iu
                        inc = 2*inc
                    }
                }
            } else {
                iu=il
                for (index in 1:n){
                    il=il-inc
                    if (il <= 1) {il=0; break}
                    else if ((xVal >= xTab[il]) == ascnd) break
                    else {
                        iu=il
                        inc = 2*inc
                    }
                }
            }
        }
            
        while(iu-il > 1){
            im=(iu+il) %/% 2
            if ((xVal >= xTab[im])== ascnd)
                il=im
            else
                iu=im
        }
        cor <- ifelse(abs(il - isav) > dj, 0,1)
        isav=il
        return(max(0,min(n-m,il-((m-2) %/% 2))))
    }


## Choose Nbr of points from Koutrouvelis tables by linear interpolation
## K for first regression : alpha and gamma
## L for second regression: beta and delta

getnbrTpoints<- function(alpha,n)
{
    xA <- as.numeric(c(1.9,1.5,1.3,1.1,0.9,0.7,0.5,0.3))
    x0 <- alpha
    if (alpha < min(xA)) x0 <- min(xA)
    if (alpha > max(xA)) x0 <- max(xA)
    yN <-  as.numeric(c(200,800,1600))
    zK <- c(9,9,10,11,11,11,22,16,14,24,18,15,28,22,18,30,24,20,86,68,56,134,124,118)
    mat <- matrix(zK,nrow=length(xA),ncol=length(yN),byrow=T)
    K <- BilinInterp_2D(xA,yN,mat,alpha,n)
    if(is.na(K))K=12
    K
}

getnbrLpoints<- function(alpha,n)
{
  xA <- as.numeric(c(1.9,1.5,1.1,0.9,0.7,0.5,0.3))
  x0 <- alpha
  if (alpha < min(xA)) x0 <- min(xA)
  if (alpha > max(xA)) x0 <- max(xA)
  yN <-  as.numeric(c(200,800,1600))
  zK <- c(9,10,11,12,14,15,16,18,17,14,14,14,24,16,16,40,38,36,70,68,66)
  mat <- matrix(zK,nrow=length(xA),ncol=length(yN),byrow=T)
  l <- BilinInterp_2D(xA,yN,mat,alpha,n)
  if(is.na(l))l=14
  l
}


