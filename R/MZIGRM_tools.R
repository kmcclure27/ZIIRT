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
#' @param a A Jx2 matrix of item slope parameters
#' @param b A JxK matrix of item intercept parameters; first column is the ZI threshold
#' @param theta An Nx2 matrix of Person Parameters
#' @param seed_person_params An optional numeric to set the seed for person parameters
#' @param seed_item_params An optional numeric to set the seed for item parameters
#' @param seed_response An optional numeric or J length numeric to set the seed for item responses
#'
#' @return A list containing the Generated Person Parameters, Item Parameters, and simulate item responses
#' @export
generate_MZIGRM_data <- function(N,J,K,rho=0,theta=NULL,a=NULL,b=NULL,
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
  #intercept0 <- matrix(stats::runif(n=J,min=-2.5,max=0),nrow=J) #b0
  #interceptK <- matrix(stats::runif(n = J*(K-1),min = -2, max = 2.5),
  #                     nrow=J)
  if(is.null(b)){
    intercept0 <- matrix(stats::runif(n=J,min=-2.5,max=0),nrow=J) #b0
    interceptK <- matrix(stats::runif(n = J*(K-1),min = -2.5, max = 2.5),
                         nrow=J)
    if(K==2){
      b = cbind(intercept0,interceptK)
    }else{
      interceptK <- t(apply(interceptK,1,sort))
    }

    b <- cbind(intercept0,interceptK)
  }else{
    stopifnot("b must have K columns for the MZI-GRM"=ncol(b)==(K),
              "b must be numeric"=is.numeric(b),
              "b must have J rows"=nrow(b)==J)
    b <- b
  }

  #if(K==2){
  #  b = cbind(intercept0,interceptK)
  #}else{
  #  interceptK <- t(apply(interceptK,1,sort))
  #  b <- cbind(intercept0,interceptK)
  #}
  item_params = list(slopes=a,intercepts=b)
  #return(item_params)
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


#' Calculates Expected Proportion of Zero Inflation
#'
#' Calculates the expected proportion of zeros due to zero inflation as a function of the presence latent variable. Assumes theta is standard bivariate normal.
#'
#' @param b0 The item intercept associated the theta0 (presence)
#' @param b1 The item intercepts associated with a zero response from theta1 (severity)
#' @param theta A vector of susceptibility/severity latent variable scores. Defaults to the population mean c(0,0)
#' @param a A vector of item discrimination parameters for susceptibility/severity latent variables. Defaults to c(1,1)
#'
#' @return The expected proportion of zeros in the MZI-GRM
#' @export

MZI_propZI <- function(b0,b1,theta=c(0,0),a=c(1,1)){
  stopifnot("b0 and b1 must be numeric"=all(is.numeric(b0),is.numeric(b1)),
            "theta must be length 2"=length(theta)==2,
            "theta must be numeric"=all(is.numeric(theta)),
            "a must be length 2"=length(a)==2,
            "a must be numeric"=all(is.numeric(a)))

  x1 = (1+exp(a[1]*theta[1]+b0))^-1
  x2 = (1+exp(a[2]*theta[2]+b1))^-1
  n = 1-x1
  d = 1-x1+x1*(1-x2)
  return(n/d)
}

#' Find b0 for MZI given a target Zero-Inflation Level
#'
#' Finds item intercept(s) for the presence latent variable in the MZI given target level of zero-inflation due to (lack of) presence and the item intercept controlling zeros from low severity
#' Assumes that the presence/severity latent variables are standard bivariate normal.
#'
#' @param propZI0 The target proportion of zeros due to (lack of) presence. Must be a decimal between 0 and 1.
#' @param b1 Item Intercept controlling potential zeros from the severity component of the MZI
#' @param n The number of candidate b0 to output. The n that provide the closest propZI0 are returned
#' @param tol A tolerance threshold that issues a warning if the resulting proportion of zero inflation is more than tol from the target
#' @param thr An alternative to n that outputs all b0 parameters with thr of the target; defaults to NULL
#'
#' @return A n-row dataframe with columns for 1) expected proportion of zero inflated due to low presence, 2) b0, and 3) b1
#' @export

MZI_find_b0 <- function(propZI0,b1,n=1,tol=0.10,thr=NULL){
  ### ERROR HANDLING HERE ###
  stopifnot("b1 must be numeric"=is.numeric(b1),
            "propZI0 must be numeric"=is.numeric(propZI0),
            "propZI0 must be positive"=propZI0>0,
            "propZI0 must be less than 1"=propZI0<1,
            "n must be an integer"=n%%1==0,
            "tol must be numeric"=is.numeric(tol),
            "tol must be less than 1"=tol<1)

  if(!is.null(thr)){
    stopifnot("thr must be numeric"=is.numeric(thr),
              "thr must be less than 1"=abs(thr)<1)
  }
  ###########################
  x = propZI_table #internal data object
  x$d = abs(propZI_table$propZI-propZI0)
  x_sub = x[x$b1==b1,] #need to handle when b0 isn't exactly in the
  if(nrow(x_sub)==0){
    y = x$b1[which.min(abs(x$b1-b1))]
    x_sub = x[x$b1==y,]
    warning(paste0("No perfect match for b1 = ",b1," found. b1 = ",y," was used"))
  }
  o = order(x_sub$d)
  x_sub = x_sub[o,];rownames(x_sub) <- NULL
  if(is.numeric(thr)){
    n = sum(x_sub$d<thr)
    message(paste0(n," possibilities within ",thr," of propZI = ",propZI0))
  }
  b0 = x_sub[1:n,c("b0","b1","propZI")]
  if(abs(x_sub$propZI[1]-propZI0)>tol){
    warning(paste0("There is a large difference between the target proportion ZI: ",propZI0,"
                    and the optimal solution: ", round(x_sub$propZI[1],4),". Proceed with caution."))
  }
  return(b0)
}

#' Find b0 for MZI given a target Zero-Inflation Level
#'
#' Finds item intercept(s) for the presence latent variable in the MZI given target level of zero-inflation due to (lack of) presence and the item intercept controlling zeros from low severity
#' Assumes that the presence/severity latent variables are standard bivariate normal.
#'
#' @param propZI0 The target proportion of zeros due to (lack of) presence. Must be a decimal between 0 and 1.
#' @param b1 Item Intercept controlling potential zeros from the severity component of the MZIGRM
#' @param n The number of candidate b0 to output. The n that provide the closest propZI0 are returned
#' @param tol A tolerance threshold that issues a warning if the resulting proportion of zero inflation is more than tol from the target
#' @param thr An alternative to n that outputs all b0 parameters with thr of the target; defaults to NULL
#'
#' @return A n-row dataframe with columns for 1) expected proportion of zero inflated due to low presence, 2) b0, and 3) b1
#' @export

MZI_find_b0 <- function(propZI0,b1,n=1,tol=0.10,thr=NULL){
  ### ERROR HANDLING  ###
  stopifnot("b1 must be numeric"=is.numeric(b1),
            "propZI0 must be numeric"=is.numeric(propZI0),
            "propZI0 must be positive"=propZI0>0,
            "propZI0 must be less than 1"=propZI0<1,
            "n must be an integer"=n%%1==0,
            "tol must be numeric"=is.numeric(tol),
            "tol must be less than 1"=tol<1)

  if(!is.null(thr)){
    stopifnot("thr must be numeric"=is.numeric(thr),
              "thr must be less than 1"=abs(thr)<1)
  }
  ###########################
  x = propZI_table #internal data object
  x$d = abs(propZI_table$propZI-propZI0)
  x_sub = x[x$b1==b1,] #need to handle when b0 isn't exactly in the
  if(nrow(x_sub)==0){
    y = x$b1[which.min(abs(x$b1-b1))]
    x_sub = x[x$b1==y,]
   # warning(paste0("No perfect match for b1 = ",b1," found. b1 = ",y," was used"))
  }
  o = order(x_sub$d)
  x_sub = x_sub[o,];rownames(x_sub) <- NULL
  if(is.numeric(thr)){
    n = sum(x_sub$d<thr)
    message(paste0(n," possibilities within ",thr," of propZI = ",propZI0))
  }
  b0 = x_sub[1:n,c("b0","b1","propZI")]
  if(abs(x_sub$propZI[1]-propZI0)>tol){
    warning(paste0("There is a large difference between the target proportion ZI: ",propZI0,"
                    and the optimal solution: ", round(x_sub$propZI[1],4),". Proceed with caution."))
  }
  return(b0)
}

#' Find b1 for MZI given a target Zero-Inflation Level
#'
#' Finds item intercept(s) for the severity latent variable in the MZI given target level of zero-inflation due to (lack of) presence and the item intercept controlling zeros from low presence
#' Assumes that the presence/severity latent variables are standard bivariate normal.
#'
#' @param propZI0 The target proportion of zeros due to (lack of) presence. Must be a decimal between 0 and 1.
#' @param b0 Item Intercept controlling potential zeros from the presence component of the MZIGRM
#' @param n The number of candidate b0 to output. The n that provide the closest propZI0 are returned
#' @param tol A tolerance threshold that issues a warning if the resulting proportion of zero inflation is more than tol from the target
#' @param thr An alternative to n that outputs all b1 parameters with thr of the target; defaults to NULL
#'
#' @return A n-row dataframe with columns for 1) expected proportion of zero inflated due to low presence, 2) b0, and 3) b1
#' @export

MZI_find_b1 <- function(propZI0,b0,n=1,tol=0.1,thr=NULL){
  ### ERROR HANDLING  ###
  stopifnot("b0 must be numeric"=is.numeric(b0),
            "propZI0 must be numeric"=is.numeric(propZI0),
            "propZI0 must be positive"=propZI0>0,
            "propZI0 must be less than 1"=propZI0<1,
            "n must be an integer"=n%%1==0,
            "tol must be numeric"=is.numeric(tol),
            "tol must be less than 1"=tol<1)

  if(!is.null(thr)){
    stopifnot("thr must be numeric"=is.numeric(thr),
              "thr must be less than 1"=abs(thr)<1)
  }
  ###########################
  x = propZI_table #internal data object
  x$d = abs(x$propZI-propZI0)
  x_sub = x[x$b0==b0,] #need to handle when b0 isn't exactly in the
  if(nrow(x_sub)==0){
    y = x$b0[which.min(abs(x$b0-b0))]
    x_sub = x[x$b0==y,]
    #warning(paste0("No perfect match for b0 = ",b0," found. b0 = ",y," was used"))
  }
  o = order(x_sub$d)
  x_sub = x_sub[o,];rownames(x_sub) <- NULL
  if(is.numeric(thr)){
    stopifnot("thr must be positive"= thr>0,
              "thr must be less than 1"= thr<1)
    n = sum(x_sub$d<thr)
    message(paste0(n," possibilities within ",thr," of propZI = ",propZI0))
  }
  b1 = x_sub[1:n,c("b0","b1","propZI")]
  #return(x_sub)
  if(abs(x_sub$propZI[1]-propZI0)>tol){
    warning(paste0("There is a large difference between the target proportion ZI: ",propZI0,"
                    and the optimal solution: ", round(x_sub$propZI[1],4),". Proceed with caution."))
  }
  return(b1)
}
