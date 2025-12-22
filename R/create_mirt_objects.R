### Attempt at a General Constructor function for an MHGRM mirt item type

#' Create a Multidimensional Hurdle Graded Response Model mirt custom item
#'
#' Creates a Multidimensional Hurdle Graded Response Model (MHGRM) with K response categories as a custom item in the mirt package
#' @param K An integer denoting the number of item response categories
#'
#' @return A mirt item from mirt::createItem() named "MHGRM_K"
#' @export
#'
create_mirt_MHGRM <- function(K){
  stopifnot("K must be an integer"=all(is.numeric(K),K))
  name=paste0("MHGRM_",K)
  par = c(a=rep(1,2),
          b=seq(-1,1,length.out=K-1))
  est = rep(TRUE,length(par))

  P.MHGRM_K <- function(par,Theta,ncat){
    D = 1.702 # Scaling Constant
    a1 = par[1]
    a2 = par[2]
    b1 = par[3]
    b = par[-c(1,2,3)]
    P0 = 1 - 1/(1+exp(-D*a1*Theta[,1]+b1))
    CP = matrix(0,nrow=nrow(Theta),ncol=length(b))
    for(k in 1:length(b)){
      CP[,k] <- 1/(1+exp(-D*a2*Theta[,2,drop=F]+b[k]))
    }
    CP = cbind(rep(1,nrow(CP)),
               CP,
               rep(0,nrow(CP)))
    P_GRM = CP[,1:(ncat-1),drop=F] - CP[,2:(ncat),drop=F]
    Pk = cbind(P0,(1-P0)*P_GRM);colnames(Pk)<-NULL
    return(Pk)
  }
  x <- mirt::createItem(
    name=name,par=par,est=est,P=P.MHGRM_K
  )
  return(x)
}

#' Create A Multidimesnional Zero-Inflated Graded Response Model mirt custom item
#'
#' NOT TESTED YET
#'
#' @param K Number of Item Response Categories (integer)
#'
#' @return A mirt custom item from mirt::createItem() named "MZIGRM_K"
#'
create_mirt_MZIGRM <- function(K){
  stopifnot("K must be an integer"=all(is.numeric(K),K))
  name=paste0("MZIGRM_",K)
  par = rep(a=c(1,2),b=seq(-1,1,length.out=K))
  est = rep(TRUE,length(par))
  P.MZIGRM_K <- function(par,Theta,ncat){
    D = 1.702 # Scaling Constant
    a1 = par[1]
    a2 = par[2]
    b0 = par[3]
    b = par[-c(1,2,3)]
    P0 = 1 - 1/(1+exp(-D*a1*Theta[,1]+b0))
    CP = matrix(0,nrow=nrow(Theta),ncol=length(b))
    for(k in 1:length(b)){
      CP[,k] <- 1/(1+exp(-D*a2*Theta[,2,drop=F]+b[k]))
    }
    CP = cbind(rep(1,nrow(CP)),
               CP,
               rep(0,nrow(CP)))
    P_GRM = CP[,1:ncat,drop=F]-CP[,2:(ncat+1),drop+F]
    Pk = cbind((P0+(1-P0)*P_GRM[,1,drop=F]),
              (1-P0)*P_GRM[,2:(ncat)+1,drop=F])
    colnames(Pk) <- NULL
    return(Pk)
  }
  x <- mirt::createItem(name=name,par=par,est=est,P=P.MZIGRM_K)
  return(x)
}
