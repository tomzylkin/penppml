#' Implementing Cross Validation
#'
#' This is the internal function called by \code{mlfitppml_int} to perform cross-validation, if the
#' option is enabled. It is available also on a stand-alone basis in case it is needed, but generally
#' users will be better served by using the wrapper \code{mlfitppml}.
#'
#' \code{xvalidate} carries out cross-validation with the user-provided IDs by holding out each one of
#' them, sequentially, as in the k-fold procedure (unless \code{testID} is specified, in which case
#' it just uses this ID for validation). After filtering out the holdout sample, the function simply
#' calls \link{penhdfeppml_int} and \link{penhdfeppml_cluster_int} to estimate the coefficients, it
#' predicts the conditional means for the held-out observations and finally it calculates the root mean
#' squared error (RMSE).
#'
#' @param testID Optional. A number indicating which ID to hold out during cross-validation. If left
#'    unspecified, the function cycles through all IDs and reports the average RMSE.
#' @param init_mu Optional: initial values of the conditional mean \eqn{\mu}, to be used as weights in the
#'     first iteration of the algorithm.
#' @param init_x Optional: initial values of the independent variables.
#' @param init_z Optional: initial values of the transformed dependent variable, to be used in the
#'     first iteration of the algorithm.
#' @param lambda Penalty parameter, to be passed on to penhdfeppml_int or penhdfeppml_cluster_int.
#' @inheritParams mlfitppml_int
#'
#' @return A list with two elements:
#' \itemize{
#'     \item \code{rmse}: root mean squared error (RMSE).
#'     \item \code{mu}: conditional means.
#' }
#' @export
#'
#' @examples
#' # First, we need to transform the data. Start by filtering the data set to keep only countries in
#' # the Americas:
#' americas <- countries$iso[countries$region == "Americas"]
#' trade <- trade[(trade$imp %in% americas) & (trade$exp %in% americas), ]
#' # Now generate the needed x, y and fes objects:
#' y <- trade$export
#' x <- data.matrix(trade[, -1:-6])
#' fes <- list(exp_time = interaction(trade$exp, trade$time),
#'             imp_time = interaction(trade$imp, trade$time),
#'             pair     = interaction(trade$exp, trade$imp))
#' # We also need to create the IDs. We split the data set by agreement, not observation:
#' id <- unique(trade[, 5])
#' nfolds <- 10
#' unique_ids <- data.frame(id = id, fold = sample(1:nfolds, size = length(id), replace = TRUE))
#' cross_ids <- merge(trade[, 5, drop = FALSE], unique_ids, by = "id", all.x = TRUE)
#' # Finally, we try xvalidate with a lasso penalty (the default) and two lambda values:
#' \dontrun{reg <- xvalidate(y = y, x = x, fes = fes, lambda = 0.001,
#'                          IDs = cross_ids$fold, verbose = TRUE)}
#'
#' @inheritSection hdfeppml_int References

xvalidate <- function(y, x, fes, IDs, testID = NULL, tol = 1e-8, hdfetol = 1e-4, colcheck_x=TRUE,colcheck_x_fes=TRUE,
                     init_mu = NULL, init_x = NULL, init_z = NULL, verbose = FALSE,
                     cluster = NULL, penalty = "lasso", method = "placeholder",
                     standardize = TRUE, penweights = rep(1, ncol(x_reg)), lambda = 0) {

  ## incorporate folds option (removed from arguments for initial release, since it wasn't doing
  # anything. Create seed option.
  if (is.null(init_mu)) {
   # init_mu <- (y + mean(y))/2
    if(is.null(fes)){
      print("fes are null")
      only_fes <- hdfeppml_int(y, fes=NULL, tol = 1e-8, hdfetol = 1e-4, colcheck_x=FALSE,colcheck_x_fes=FALSE, mu = NULL, saveX = TRUE,
                             init_z = NULL, verbose = FALSE, maxiter = 1000, cluster = NULL, vcv = TRUE)
      init_mu <- only_fes$mu
    } else {
      only_fes <- hdfeppml_int(y, fes=fes, tol = 1e-8, hdfetol = 1e-4, colcheck_x=FALSE,colcheck_x_fes=FALSE, mu = NULL, saveX = TRUE,
                               init_z = NULL, verbose = FALSE, maxiter = 1000, cluster = NULL, vcv = TRUE)
      mu <- only_fes$mu
    }
  }
  if (is.null(init_x)) {
    init_x <- x
  }
  x <- data.matrix(x)
  init_x <- data.matrix(init_x)

  uniq_IDs <- levels(factor(IDs))
  n_IDs    <- length(uniq_IDs)
  mu <- init_mu
  # < Initialize mu, x

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

  # print("length ids and omitid")
  # print(length(IDs))
  #drop 1 id at a time and predict its mean values out of sample
  for (i in start_loop:n_IDs) {
#    print(i)
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

    # this is all to create some temporary fixed effects and fixed effects names
    if(is.null(fes)){
      fes_temp <- NULL
    }else{
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
    }

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
    if (method == "plugin"){
      cluster_reg    <- cluster[insample]
      plugin_xval <- penhdfeppml_cluster_int(y=y_temp,x=x_reg,fes=fes_temp,lambda=lambda,cluster=cluster_reg,tol=tol,hdfetol=hdfetol,verbose=FALSE, colcheck_x=FALSE,colcheck_x_fes=FALSE)
      if(verbose) {
        print("hdfeppml finished")
      }
      mu_temp <- plugin_xval$mu
      b  <-plugin_xval$beta
    }
    else if (penalty == "ridge"){
      ridge_xval <-  penhdfeppml_int(y=y_temp,x=x_reg,fes=fes_temp,lambda=lambda,tol=tol,hdfetol=hdfetol,penalty="ridge",standardize=standardize,verbose=verbose,colcheck_x=FALSE,colcheck_x_fes=FALSE)
      if(verbose) {
        print("hdfeppml finished")
      }
      mu_temp <- ridge_xval$mu
      b  <-ridge_xval$beta

    } else {
      lasso_xval <- penhdfeppml_int(y=y_temp,x=x_reg,fes=fes_temp,lambda=lambda,tol=tol,hdfetol=hdfetol,penalty="lasso",
                                standardize=standardize,verbose=verbose,penweights=penweights,colcheck_x=FALSE,colcheck_x_fes=FALSE)

      if(verbose) {
       print("hdfeppml finished")
      }
      mu_temp <- lasso_xval$mu
      b  <-lasso_xval$beta
    }
    # < Calculate coefficients b with subsample

    mu_temp <- compute_fes(y=y,fes=fes,x=x,b=b,insample_obs=(IDs!=omitID),onlymus=TRUE,tol=tol,verbose=verbose)
    mu[which(IDs==omitID)] <- mu_temp[which(IDs==omitID)]

    if (verbose) {
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

#' Filtering fixed effect lists
#'
#' A helper function for \code{xvalidate} that filters a list of fixed effects and returns the modified
#' list. Used to split the fixed effects for cross-validation.
#'
#' @param fe_list A list of fixed effects.
#' @param select_obs A vector of selected observations / rows.
#' @param list Logical. If \code{TRUE}, it returns a list. Otherwise, a data frame.
#'
#' @return A modified list of fixed effects.


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
      else {
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
#' This function is a helper for \code{xvalidate} that computes FEs using PPML First Order Conditions
#' (FOCs).
#'
#' @param y Dependent variable (a vector).
#' @param fes List of fixed effects.
#' @param x Regressor matrix.
#' @param b A vector of coefficient estimates.
#' @param insample_obs Vector of observations used to estimate the \code{b} coefficients..
#' @param onlymus Logical. If \code{TRUE}, returns only the conditional means.
#' @param tol A tolerance parameter.
#' @param verbose Logical. If \code{TRUE}, prints messages to the console while evaluating.
#'
#' @return If \code{onlymus = TRUE}, the vector of conditional means. Otherwise, a list with two
#' elements:
#' \itemize{
#'     \item \code{mu}: conditional means.
#'     \item \code{fe_values}: fixed effects.
#' }

compute_fes <- function(y, fes, x, b, insample_obs = rep(1, n),
                        onlymus = FALSE, tol = 1e-8, verbose = FALSE) {

  n <- length(y)

  # recover FEs
  # Initialize critical values and loop value
  crit <- 1
  j <- 0
  tol <- tol
  deviance <-  0

  # If coefficient not selected set to zero
  b[which(is.na(b))] <- 0
  if(is.null(fes)){
    mus <- data.frame(intercept=1, mu = exp(x %*% b) * (insample_obs), y = y * insample_obs)
  } else {
    mus <- data.frame(fes, mu = exp(x %*% b) * (insample_obs), y = y * insample_obs)
  }

  if(verbose) {
    print("recovering FEs")
  }

  # iteratively solve for FEs
  if(is.null(fes)){
    while (crit>tol) {
        if (j==0) {
          intercept <- rep(1,n) #initialize (exponentiated) FEs using 1's
         # names(mus)["intercept"] <- 1        #change name of FE id within mus to a generic id
        }
        mus <- within(mus, {y_sum  = sum(y)} )
        mus <- within(mus, {mu_sum = sum(mu)} )
        print(head(mus))
        # update FEs using PPML FOCs
        adj <- (mus$y_sum / mus$mu_sum)
        adj[which(mus$y_sum==0)]<-0
        if(onlymus==TRUE){
          intercept <- intercept * adj
        }
        mus$mu <- mus$mu * adj

      # compute deviance and verify convergence
      last_dev <- deviance
      temp <- mus$y * log(mus$y/mus$mu) - (mus$y-mus$mu)
      temp[which(mus$y==0)] <- -mus$mu[which(mus$y==0)]
      deviance <- 2 * sum(temp)
      denom_eps = max( min(deviance, last_dev) , 0.1 )
      crit <- abs(deviance-last_dev)/denom_eps

      j <- j + 1
    } # < end while
  } else {
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
  }

  # return values
  mu <- exp(x%*%b)
  if(!is.null(fes)){
  for (f in 1:length(fes)) {
    fe_value <- paste("fe_value",f,sep="")
    mu <- mu * get(fe_value)
  }
  }else{
    mu <- mu * intercept
  }
  if (onlymus) {
    return(mu)
  }
  else{
    if(!is.null(fes)){
    fe_values <- matrix(nrow=n,ncol=length(fes))
    fe_values <- data.frame(fe_values)
    for (f in 1:length(fes)) {
      fe_name  <- paste("fe",f,sep="")
      fe_value <- paste("fe_value",f,sep="")
      fe_values[,f]    <- get(fe_value)
      names(fe_values)[f] <- fe_name
    }
    returnlist = list("mu"=mu,"fe_values"=fe_values)
    } else {
    returnlist = list("mu"=mu,"fe_values"=intercept)
    }
  }
}
