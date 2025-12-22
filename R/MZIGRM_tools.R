

#' MZIGRM Response Probabilities
#'
#' Calculates the MZI Graded Response Model response category probabilities
#'
#' @param theta An N x 2 Matrix of Person Parameters
#' @param slope A 1 x 2 Matrix of item slopes
#' @param intercept A 1 x K Matrix of Item Intercepts. The first intercept corresponds to the 2PL component
#' @param D A scaling Constant
#'
#' @return An N x K matrix of response probabilities
#'
MZIGRM_mat2 <- function(theta,slope,intercept,D=1.702){
  if(is.vector(theta)){
    theta=matrix(theta,nrow=1)
  }
  N = nrow(theta)
  n_ints = ncol(intercept)
  P0 = GRM_mat2(theta=theta[,1,drop=F],
                slope=slope[,1,drop=F],
                intercept=intercept[,1,drop=F],D=D)
  P_GRM = GRM_mat2(theta=theta[,2,drop=F],
                   slope=slope[,2,drop=F],
                   intercept=intercept[,2:n_ints,drop=F],D=D)
  Pk = cbind(P0[,1]+(P0[,2])*P_GRM[,1],
             (P0[,2]*P_GRM[,2:n_ints]))
  return(Pk)
}

#' Simulate Item Responses for the MZIGRM
#'
#' Simulates Item Responses for the multidimensional zero-inflated graded response model. This function wraps MHGRM_simulate_response
#'
#' @inheritParams MHGRM_simulate_response
#'
#' @return An N length vector of item responses 0:(K-1)
#'
MZIGRM_simulate_response <- function(x,seed=NULL){
  y = MHGRM_simulate_response(x=x,seed=seed)
  return(y)
}

#' Generate MZIGRM Data
#'
#' Generates Some Factory MZIGRM Item Parameters, Person Parameters, and Item Responses
#'
#' @param N Number of Observations
#' @param J Number of Items
#' @param K Number of Response Categories. Must be greater than 1 for the MZIGRM
#' @param rho Correlation Between Person Parameters
#' @param seed_person_params An optional numeric to set the seed for person parameters
#' @param seed_item_params An optional numeric to set the seed for item parameters
#' @param seed_response An optional numeric or J length numeric to set the seed for item responses
#'
#' @return A list containing the Generated Person Parameters, Item Parameters, and simulate item responses
#' @export
generate_MZIGRM_data <- function(N,J,K,rho=0,
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
  pps <- MASS::mvrnorm(n=N,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),nrow=2,byrow=T))
  if(!is.null(seed_item_params)){
    stopifnot("seed_item_params should be a scalar"=length(seed_item_params)==1,
              "seed_item_params must be numeric"=is.numeric(seed_item_params))
    set.seed(seed_item_params)
  }
  a <- matrix(stats::runif(n=2*J,min=1,max=4),nrow=J,ncol=2)
  intercept0 <- matrix(stats::runif(n=J,min=-2.5,max=0),nrow=J) #b0
  interceptK <- matrix(stats::runif(n = J*(K-1),min = -2.5, max = 2.5),
                       nrow=J)
  if(K==2){
    b = cbind(intercept0,interceptK)
  }else{
    interceptK <- t(apply(interceptK,1,sort))
    b <- cbind(intercept0,interceptK)
  }
  item_params = list(slopes=a,intercepts=b)

  resp_ps <- array(data=NA, dim = c(N,K,J))
  for(j in 1:J){
    resp_ps[,,j] <- MZIGRM_mat2(theta=pps,
                               slope=a[j,,drop=F],
                               intercept=b[j,,drop=F])
  }
  item_resps <- matrix(NA,nrow=N,ncol=J)
  if(!is.null(seed_response)){
    if(length(seed_response)==1){
      stopifnot("seed_response must be numeric"=is.numeric(seed_response))
      set.seed(seed_response)
      for(j in 1:J){
        item_resps[,j] <- MZIGRM_simulate_response(x=resp_ps[,,j])
      }
    }else if(length(seed_response)==J){
      stopifnot("seed_response must be numeric"=is.numeric(seed_response))
      for(j in 1:J){
        item_resps[,j] <- MZIGRM_simulate_response(x=resp_ps[,,j],
                                                  seed=seed_response[j])
      }
    }else{
      warning("seed_response must be a scalar or J length numeric vector. seed not set")
      for(j in 1:J){
        item_resps[,j] <- MZIGRM_simulate_response(x=resp_ps[,,j],
                                                  seed=seed_response[j])
      }
    }
  }else{
    for(j in 1:J){
      item_resps[,j] <- MZIGRM_simulate_response(x=resp_ps[,,j])
    }
  }
  colnames(item_resps) <- paste0("x",1:J)
  out = list(item_params=item_params,
             person_params = pps,
             item_resps=item_resps)
  return(out)
}
