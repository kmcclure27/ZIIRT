#### Multidimensional Hurdle Graded Response Model Functions ####

#' GRM Item Response Probabilities
#'
#' Calculates the Graded Response Model response category probabilities
#'
#' @param theta An N x M matrix of person parameters
#' @param slope A 1 x M matrix of item slopes
#' @param intercept A 1 x K-1 matrix of item intercepts
#' @param D A scaling Constant
#'
#' @return An N x K matrix of item response probabilities
#'
GRM_mat2 <- function(theta,slope,intercept,D=1.702){
  if(is.vector(theta)){
    theta=matrix(theta,nrow=1)
  }
  N = nrow(theta)
  n_ints = length(intercept)
  eta <- matrix(-D*theta%*%t(slope),nrow=N,
                ncol=n_ints)
  int_mat = matrix(as.numeric(intercept),
                   nrow=N,ncol=n_ints,
                   byrow=T)
  y = eta + int_mat
  cum_prob = cbind(rep(1,N),
                   (1+exp(y))^-1,
                   rep(0,N))
  p = cum_prob[,1:(n_ints+1)] - cum_prob[,2:(n_ints+2)]
  return(p)
}

#' MHGRM Response Probabilities
#'
#' Calculates the MH Graded Response Model response category probabilities
#'
#' @param theta An N x 2 Matrix of Person Parameters
#' @param slope A 1 x 2 Matrix of item slopes
#' @param intercept A 1 x K-1 Matrix of Item Intercepts
#' @param D A scaling Constant
#'
#' @return An N x K matrix of response probabilities
#'
MHGRM_mat2 <- function(theta,slope,intercept,D=1.702){
  if(is.vector(theta)){
    theta=matrix(theta,nrow=1)
  }
  N = nrow(theta)
  n_ints = ncol(intercept)
  P0 = GRM_mat2(theta=theta[,1,drop=F],
                slope=slope[,1,drop=F],
                intercept=intercept[,1,drop=F],D=D)
  if(ncol(intercept)==1){
    Pk=P0
    return(Pk)
  }
  P_GRM = GRM_mat2(theta=theta[,2,drop=F],
                   slope=slope[,2,drop=F],
                   intercept=intercept[,2:n_ints,drop=F],D=D)
  Pk = cbind(P0[,1],
             P_GRM*P0[,2])
  return(Pk)
}

#' Simulate Item Responses from MHGRM
#'
#' Simulates item responses under the MHGRM give an N x K matrix of response probabilities
#'
#' @param x An N x K matrix of response probabilities
#' @param seed A numeric value to set the random seed
#'
#' @return An N length vector of item responses ranging from 0:(K-1)
#'
MHGRM_simulate_response <- function(x,seed=NULL){
  if(!is.null(seed)){
    stopifnot("seed must be numeric"=is.numeric(seed))
    set.seed(seed)
  }
  y = apply(x,1,function(x){
    sample(1:length(x),size=1,prob=x)-1
  })
  return(y)
}


#' Generate MHGRM Data
#'
#' Generates Some Factory MHGRM Item Parameters, Person Parameters, and Item Responses
#'
#' @param N Number of Observations
#' @param J Number of Items
#' @param K Number of Response Categories. Must be greater than 2 for the MHGRM
#' @param rho Correlation Between Person Parameters
#' @param theta An Nx2 matrix of Person Parameters
#' @param a A Jx2 matrix of slope parameters
#' @param b A JxK-1 matrix of intercept parameters
#' @param seed_person_params An optional numeric to set the seed for person parameters
#' @param seed_item_params An optional numeric to set the seed for item parameters
#' @param seed_response An optional numeric or J length numeric to set the seed for item responses
#'
#' @return A list containing the Generated Person Parameters, Item Parameters, and simulate item responses
#' @export
#'
generate_MHGRM_data <- function(N,J,K,rho=0,theta=NULL,a=NULL,b=NULL,
                                seed_person_params=NULL,seed_item_params=NULL,
                                seed_response=NULL){
  stopifnot("K must be greater than 1" = K>1,
            "K must be an integer" = K%%1==0,
            "N must be an integer" = N%%1==0,
            "J must be an integer" = J%%1==0,
            "rho must be numeric between -1 and 1"=all(is.numeric(rho),
                                                       length(rho)==1,
                                                       abs(rho)<1))


  if(!is.null(seed_person_params)){
    stopifnot("seed_person_params should be a scalar"=length(seed_person_params)==1,
              "seed_person_params must be numeric"=is.numeric(seed_person_params))
    set.seed(seed_person_params)
  }
  if(is.null(theta)){
    pps <- MASS::mvrnorm(n=N,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),nrow=2,byrow=T))
  }else{
    stopifnot("theta must have two columns" = ncol(theta)==2,
              "theta must have N rows"=nrow(theta)==N,
              "theta must be numeric"=is.numeric(theta))
    pps<-theta
  }
  if(!is.null(seed_item_params)){
    stopifnot("seed_item_params should be a scalar"=length(seed_item_params)==1,
              "seed_item_params must be numeric"=is.numeric(seed_item_params))
    set.seed(seed_item_params)
  }
  if(is.null(a)){
    a <- matrix(stats::runif(n=2*J,min=1,max=4),nrow=J,ncol=2)
  }else{
    stopifnot("a must have two columns"=ncol(a)==2,
              "a must be numeric"=is.numeric(a),
              "a must have J rows"=nrow(a)==J)
    a <- a
  }
  if(is.null(b)){
    intercept0 <- matrix(stats::runif(n=J,min=-2.5,max=0),nrow=J) #b0
    interceptK <- matrix(stats::runif(n = J*(K-2),min = -2.5, max = 2.5),
                         nrow=J)

    if(K==3){
      b = cbind(intercept0,interceptK)
    }else{
      interceptK <- t(apply(interceptK,1,sort))
    }
    if(K==2){
      b = matrix(intercept0,ncol=1)
    }else{
      b <- cbind(intercept0,interceptK)
    }
  }else{
    stopifnot("b must have K-1 columns for the MH-GRM"=ncol(b)==(K-1),
              "b must be numeric"=is.numeric(b),
              "b must have J rows"=nrow(b)==J)
    b <- b
  }

  item_params = list(slopes=a,intercepts=b)
  resp_ps <- array(data=NA, dim = c(N,K,J))
  for(j in 1:J){
    resp_ps[,,j] <- MHGRM_mat2(theta=pps,
                               slope=a[j,,drop=F],
                               intercept=b[j,,drop=F])
  }
  item_resps <- matrix(NA,nrow=N,ncol=J)
  if(!is.null(seed_response)){
    if(length(seed_response)==1){
      stopifnot("seed_response must be numeric"=is.numeric(seed_response))
      set.seed(seed_response)
      for(j in 1:J){
        item_resps[,j] <- MHGRM_simulate_response(x=resp_ps[,,j])
      }
    }else if(length(seed_response)==J){
      stopifnot("seed_response must be numeric"=is.numeric(seed_response))
      for(j in 1:J){
        item_resps[,j] <- MHGRM_simulate_response(x=resp_ps[,,j],
                                                  seed=seed_response[j])
      }
    }else{
      warning("seed_response must be a scalar or J length numeric vector. seed not set")
      for(j in 1:J){
        item_resps[,j] <- MHGRM_simulate_response(x=resp_ps[,,j],
                                                  seed=seed_response[j])
      }
    }
  }else{
    for(j in 1:J){
      item_resps[,j] <- MHGRM_simulate_response(x=resp_ps[,,j])
    }
  }
  colnames(item_resps) <- paste0("x",1:J)
  out = list(item_params=item_params,
             person_params = pps,
             item_resps=item_resps)
  return(out)
}
