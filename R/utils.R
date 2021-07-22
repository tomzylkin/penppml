## !!! CURRENTLY NOT USING !!!!
# setupLambda <- function(X, y, lambda.min, nlambda, penalty.factor) {
#  n <- nrow(X)
#  p <- ncol(X)

#  ## Determine lambda.max
#  ind <- which(penalty.factor!=0)
#  if (length(ind)!=p) {
#    fit <- glm(y~X[, -ind], family="gaussian")
#  } else {
#    fit <- glm(y~1, family=family)
#  }

#  zmax <- .Call("maxprod", X, fit$residuals, ind, penalty.factor) / n
#  lambda.max <- zmax/alpha

#  if (lambda.min==0) {
#    lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)), 0)
#  } else {
#    lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), len=nlambda))
#  }
#
#  if (length(ind)!=p) lambda[1] <- lambda[1] * 1.000001
#  lambda
#}


#' Title
#'
#' @param y
#' @param x
#' @param fes
#' @param hdfetol
#' @param selectobs
#'
#' @return
#' @export
#'
#' @examples

collinearity_check <- function(y,x,fes,hdfetol,selectobs=NULL) {
  # selectobs doesn't currently do anything.

  mu  <- (y + mean(y))/2
  z   <- (y-mu)/mu + log(mu)
  reg_z  <- matrix(z)
  reg_x  <- x
  mu  <- (y + mean(y))/2
  z_resid <- demeanlist(reg_z,fes,weights=sqrt(mu),eps=hdfetol)
  x_resid <- demeanlist(reg_x,fes,weights=sqrt(mu),eps=hdfetol)

  check <- lm.wfit(x_resid, z_resid, mu) #faster than lmfit
  check$coefficients
  include_x <- which(!is.na(check$coefficients))
}


#' Title
#'
#' @param e
#' @param cluster
#' @param x
#'
#' @return
#' @export
#'
#' @examples

cluster_matrix <- function(e,cluster,x) {
  K <- ncol(x)
  vars      <- data.frame(e=e,cluster=factor(cluster,exclude=TRUE),x=x)
  vars      <- vars[order(vars$cluster),] #5780 is missing...
  vars$indx <- with(vars, ave(seq_along(cluster), cluster, FUN=seq_along))
  vars      <- complete(vars, cluster, indx, fill = list(e=0,x=0))
  vars      <- vars %>% fill(everything(),0)

  T <- max(vars$indx)

  if(is.null(K)){
    X <- data.matrix(vars[,4])
  } else{
    X <- data.matrix(vars[,4:(3+K)])
  }

  e <- as.matrix(vars[,3])
  rm(vars)

  ee   <- manyouter(e,e,T)  #many outer products for vectors of length T
  XeeX <- xeex(X,ee)        #XeeX matrix for cluster-robust SEs
  return(XeeX)
}


#' Title
#'
#' @param x
#' @param weights
#' @param intercept
#' @param return.sd
#'
#' @return
#' @export
#'
#' @examples

standardize_wt <- function(x,weights=rep(1/n,n),intercept=TRUE,return.sd=FALSE){
  n     <- nrow(x)
  nvars <- ncol(x)
  if (intercept) {
    xm <- fastwmean(x,weights)
  } else {
    xm <- rep(0.0, times = nvars)
  }
  xs <- faststddev(x,weights)

  if (return.sd==TRUE) {
    return(xs)
  }
  else{
    if (!inherits(x, "sparseMatrix")) x_std <- t((t(x) - xm) / xs)
    return(x_std)
  }
}


#' Title
#'
#' @param x
#' @param y
#' @param weights
#' @param lambda
#' @param standardize
#'
#' @return
#' @export
#'
#' @examples

fastridge <- function(x,y,weights=rep(1/n,n),lambda,standardize=TRUE){
  n <- length(y)
  if(standardize){
    b <- fastridgeCpp(sqrt(weights)*standardize_wt(x,weights),sqrt(weights)*y,lambda)
    beta <- b * faststddev(x,weights)
  } else{
    beta <- fastridgeCpp(sqrt(weights)*x,sqrt(weights)*y,lambda)
  }
  result <- list("beta" = beta)
  return(result)
}
