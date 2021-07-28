#' Implementing Cross Validation
#'
#' TODO: put some description here.
#'
#' @param testID TODO: check this.
#' @param init_mu TODO: check this.
#' @param init_x TODO: check this.
#' @param init_z TODO: check this.
#' @param folds TODO: check this.
#' @param lambda Penalty parameter, to be passed on to penhdfeppml or penhdfeppml_cluster.
#' @inheritParams mlfitppml
#'
#' @return A list (TODO: check this).
#' @export
#'
#' @examples #TODO: add examples.

xvalidate = function(y, x, fes, IDs, testID = NULL, tol = 1e-8, hdfetol = 1e-4, colcheck = TRUE,
                     init_mu = NULL, init_x = NULL, init_z = NULL, verbose = FALSE, folds = NULL,
                     cluster = NULL, penalty = c("none", "lasso", "ridge"), method = "placeholder",
                     standardize = TRUE, penweights = rep(1, ncol(x_reg)), lambda = 0) {

  set.seed(1)

  ## incorporate folds option. Create seed option.
  if (is.null(init_mu)) {
    init_mu <- (y + mean(y))/2
  }
  if (is.null(init_x)) {
    init_x <- x
  }
  x <- data.matrix(x)
  init_x <- data.matrix(init_x)

  uniq_IDs <- levels(factor(IDs))
  n_IDs    <- length(uniq_IDs)
  mu <- init_mu

  #ID = 0 assumed to be omitted group
  if (min(IDs)==0) {
    start_loop <- 2
  }
  else {
    start_loop <- 1
  }
  if(!is.null(testID)) {
    start_loop <- testID
    n_IDs      <- testID
  }

  #drop 1 id at a time and predict its mean values out of sample
  for (i in start_loop:n_IDs) {

    omitID <- uniq_IDs[i]
    if(!is.null(testID)) {
      omitID <- testID
    }
    if (verbose==TRUE) {
      print(i)
      print(omitID)
    }
    insample <- which(IDs!=omitID)
    if(verbose){
      print(length(insample))
      print(length(y))
    }

    #select y and x
    y_temp   <- y[insample]
    x_temp   <- x[insample,]
    x_reg    <- init_x[insample,]

    # this is all to create some temporary fixed efects and fixed effects names
    fes_temp <- select_fes(fes,insample,list=FALSE)
    for (f in 1:length(fes)) {
      fe_name <- paste("fe",f,sep="")
      assign(fe_name,factor(fes_temp[,f]))
      if (f == 1) {
        temp <- list(fe1=factor(fes_temp[,f]))
      }
      else{
        temp[[paste("fe",f,sep="")]] <- factor(fes_temp[,f])
      }
    }
    fes_temp <- temp
    rm(temp)


    # initialize mu and z
    mu_temp  <- init_mu[insample]
    z_temp   <- (y_temp-mu_temp)/mu_temp + log(mu_temp)

    if (is.null(init_z)) {
      z_reg <- z_temp
    }
    else{
      z_reg <- init_z[insample]
    }

    if(verbose) {
      print("beginning hdfeppml")
    }

    # may not need colcheck? really just need mu here.
    if (method == "iterative"){
      cluster_reg    <- cluster[insample]
      plugin_xval <- penhdfeppml_cluster(y=y_temp,x=x_reg,fes=fes_temp,lambda=lambda,cluster=cluster_reg,tol=tol,hdfetol=hdfetol,verbose=FALSE)
      if(verbose) {
        print("hdfeppml finished")
      }
      mu_temp <- plugin_xval$mu
      b  <-plugin_xval$beta
    }
    else if(penalty=="ridge"){
      ridge_xval <-  penhdfeppml(y=y_temp,x=x_reg,fes=fes_temp,lambda=lambda,tol=tol,hdfetol=hdfetol,penalty="ridge",standardize=standardize,verbose=verbose)
      if(verbose) {
        print("hdfeppml finished")
      }
      mu_temp <- ridge_xval$mu
      b  <-ridge_xval$beta

    } else {
      lasso_xval <- penhdfeppml(y=y_temp,x=x_reg,fes=fes_temp,lambda=lambda,tol=tol,hdfetol=hdfetol,penalty="lasso",
                                standardize=standardize,verbose=verbose,penweights=penweights)

      if(verbose) {
       print("hdfeppml finished")
      }
      mu_temp <- lasso_xval$mu
      b  <-lasso_xval$beta
    }

    mu_temp <- compute_fes(y=y,fes=fes,x=x,b=b,insample_obs=(IDs!=omitID),onlymus=TRUE,tol=tol,verbose=verbose)
    mu[which(IDs==omitID)] <- mu_temp[which(IDs==omitID)]

    if(verbose) {
      print("assigned FEs and computed means")
    }
  }

  y_temp  <- y[which(IDs!=0)]
  mu_temp <- mu[which(IDs!=0)]

  # rmse (for FTA pairs)
  ss <- sum((y_temp-mu_temp)^2)
  rmse <- sqrt(sum(ss)/length(ss))/stats::sd(y)

  # calculate deviance
  #temp <-  -(y_temp * log(y_temp/mu_temp) - (y_temp-mu_temp))  #mu_temp can be zero without y_temp being zero, in which case deviance is undefined.
  #temp[which(y_temp==0)] <- -mu[which(y_temp==0)]
  #deviance <- -2 * sum(temp)/length(y_temp)

  if(verbose){
    print(rmse)
    #print(min(temp))
    #print(y[which(temp==-Inf)])
    #print(deviance)
  }
  #test_output <- cbind(y_temp,mu_temp,temp)

  returnlist <- list("rmse"=rmse,"mu"=mu)
}

#' Title
#'
#' A helper function for \code{xvalidate} (TODO: add some description of what this does).
#'
#' @param fe_list A list of fixed effects.
#' @param select_obs A vector of selected observations / rows.
#' @param list Logical (TODO: check what this does).
#'
#' @return A modified list of fixed effects.
#'
#' @examples #TODO: add some examples.

select_fes <- function(fe_list, select_obs, list = TRUE) {
  fes_temp <- data.frame(fe_list)
  fes_temp <- fes_temp[select_obs,]
  if (list) {
    for (f in 1:length(fe_list)) {
      fe_name <- paste("fe", f, sep = "")
      assign(fe_name, factor(fes_temp[, f]))
      if (f == 1) {
        temp <- list(fe1 = factor(fes_temp[, f]))
      }
      else{
        temp[[paste("fe", f, sep = "")]] <- factor(fes_temp[, f])
      }
    }
    fes_temp <- temp
    rm(temp)
  }
  return(fes_temp) #new fes in list form
}


#' Fixed Effects Computation
#'
#' Compute FEs using PPML FOCs (TODO: explain this better?).
#'
#' @param y Dependent variable (a vector).
#' @param fes List of fixed effects.
#' @param x Regressor matrix.
#' @param b TODO: check what this does.
#' @param insample_obs Numeric vector with selected observations (TODO: check this better).
#' @param onlymus Logical. If \code{TRUE}, returns only mus (TODO: explain this).
#' @param tol A tolerance parameter.
#' @param verbose Logical. If \code{TRUE}, prints messages to the console while evaluating.
#'
#' @return A list (TODO: check and explain better).
#'
#' @examples #TODO: add some examples.

compute_fes <- function(y, fes, x, b, insample_obs = rep(1, n),
                        onlymus = FALSE, tol = 1e-8, verbose = FALSE) {

  n <- length(y)

  # recover FEs
  crit <- 1
  j <- 0
  tol <- tol
  deviance = 0

  b[which(is.na(b))]<-0
  mus <- data.frame(fes,mu=exp(x%*%b)*(insample_obs),y=y*insample_obs)

  if(verbose) {
    print("recovering FEs")
  }

  # iteratively solve for FEs
  while (crit>tol) {
    for (f in 1:length(fes)) {
      fe_name  <- paste("fe",f,sep="")
      fe_value <- paste("fe_value",f,sep="")
      if (j==0) {
        assign(fe_value,rep(1,n)) #initialize (exponentiated) FEs using 1's
        names(mus)[f] <- fe_name          #change name of FE id within mus to a generic id
      }
      mus <- within(mus, {y_sum  = stats::ave(y,get(fe_name),FUN=sum)} )
      mus <- within(mus, {mu_sum = stats::ave(mu,get(fe_name),FUN=sum)} )

      # update FEs using PPML FOCs
      adj <- (mus$y_sum / mus$mu_sum)
      adj[which(mus$y_sum==0)]<-0
      assign(fe_value,get(fe_value)*adj)
      mus$mu <- mus$mu * adj
    }

    # compute deviance and verify convergence
    last_dev <- deviance
    temp <- mus$y * log(mus$y/mus$mu) - (mus$y-mus$mu)
    temp[which(mus$y==0)] <- -mus$mu[which(mus$y==0)]
    deviance <- 2 * sum(temp)
    denom_eps = max( min(deviance, last_dev) , 0.1 )
    crit <- abs(deviance-last_dev)/denom_eps

    j <- j + 1
  }

  # return values
  mu <- exp(x%*%b)
  for (f in 1:length(fes)) {
    fe_value <- paste("fe_value",f,sep="")
    mu <- mu * get(fe_value)
  }
  if(onlymus) {
    return(mu)
  }
  else{
    fe_values <- matrix(nrow=n,ncol=length(fes))
    fe_values <- data.frame(fe_values)
    for (f in 1:length(fes)) {
      fe_name  <- paste("fe",f,sep="")
      fe_value <- paste("fe_value",f,sep="")
      fe_values[,f]    <- get(fe_value)
      names(fe_values)[f] <- fe_name
    }
    returnlist = list("mu"=mu,"fe_values"=fe_values)
  }
}
