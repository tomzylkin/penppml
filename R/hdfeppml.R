#' Title
#'
#' @param y
#' @param x
#' @param fes
#' @param tol
#' @param hdfetol
#' @param colcheck
#' @param selectobs
#' @param mu
#' @param saveX
#' @param init_z
#' @param verbose
#' @param maxiter
#' @param cluster
#' @param vcv
#'
#' @return
#' @export
#'
#' @examples

hdfeppml <- function(y, x, fes, tol = 1e-8, hdfetol = 1e-4, colcheck = TRUE, selectobs = NULL,
                     mu = NULL, saveX = TRUE, init_z = NULL, verbose = FALSE, maxiter = 1000,
                     cluster = NULL, vcv = TRUE) {

  if (is.null(selectobs)) {
    selectobs <- rep(TRUE,length(y))
  }
  else {
    #y <- y[selectobs]
    #x <- x[selectobs]
  }

  # number of observations (needed for deviance)
  n   <- length(selectobs)

  x<-data.matrix(x)

  # estimation algorithm
  crit <- 1
  iter <-0
  old_deviance <- 0
  include_x <- 1:ncol(x)

  b <- matrix(NA,nrow=ncol(x),ncol=1)
  xnames <- colnames(x)
  if (colcheck==TRUE){
    if (verbose == TRUE) {
      print("checking collinearity")
    }
    include_x <- collinearity_check(y,x,fes,1e-6)
    x <- x[,include_x]
  }

  if (verbose == TRUE) {
    print("beginning estimation")
  }

  while (crit>tol & iter<maxiter) {
    iter<-iter+1

    if (verbose == TRUE) {
      print(iter)
    }

    if(iter==1) {
      ## initilize "mu"
      y   <- y[selectobs]
      if (is.null(mu)) mu  <- (y + mean(y))/2
      z   <- (y-mu)/mu + log(mu)
      eta <- log(mu)
      last_z <- z
      if (is.null(init_z)) {
        reg_z  <- matrix(z)
      } else{
        reg_z <- init_z
      }
      #reg_x  <- x[selectobs,]  ## this won't work if x has one column. Will work otherwise.
      reg_x  <- x

      ## how to subset fes?
    } else {
      last_z <- z
      z <- (y-mu)/mu + log(mu)
      reg_z  <- matrix(z - last_z + z_resid)
      reg_x  <- x_resid
      ## colnames(reg_x)   <- colnames(x)
    }
    if (verbose == TRUE) {
      print("within transformation step")
    }
    z_resid <- demeanlist(reg_z,fes,weights=sqrt(mu),eps=hdfetol)
    x_resid <- demeanlist(reg_x,fes,weights=sqrt(mu),eps=hdfetol)

    if (verbose == TRUE) {
      print("obtaining coefficients")
    }
    reg <- fastolsCpp(sqrt(mu)*x_resid,sqrt(mu)*z_resid)  #faster without constant
    b[include_x] <-  reg
    reg <- list("coefficients"=b) # this rewrites reg each time. Better to initialize upfront?

    if (verbose == TRUE) {
      print(iter)
      #print(b)
    }

    if (length(include_x)==1) {
      reg$residuals <- z_resid - x_resid * b[include_x]
    } else {
      reg$residuals <- z_resid - x_resid %*% b[include_x]
    }
    mu <- as.numeric(exp(z[selectobs]-reg$residuals))
    if (verbose == TRUE) {
      print("info on residuals")
      print(max(reg$residuals))
      print(min(reg$residuals))

      print("info on means")
      print(max(mu))
      print(min(mu))

      print("info on coefficients")
      print(max(b[include_x]))
      print(min(b[include_x]))
    }

    if (verbose == TRUE) {
      print("calculating deviance")
    }

    # calculate deviance
    temp <-  -(y * log(y/mu) - (y-mu))
    temp[which(y==0)] <- -mu[which(y==0)]
    if (!missing(selectobs)){
      temp[which(!selectobs)] <-0
    }

    deviance <- -2 * sum(temp)/n

    if(deviance<0) deviance = 0

    delta_deviance <- old_deviance - deviance

    if (!is.na(delta_deviance) & (deviance < 0.1 * delta_deviance)) {
      delta_deviance = deviance
    }
    if (verbose == TRUE) {
      print("checking critical value")
    }
    denom_crit = max(c( min(c(deviance, old_deviance)) , 0.1 ))
    crit = abs(delta_deviance) / denom_crit
    if (verbose == TRUE) {
      print(deviance)
      print(crit)
    }
    old_deviance <- deviance
  }

  temp <-  -(y * log(y/mu) - (y-mu))
  temp[which(y==0)] <- 0
  if (!missing(selectobs)){
    temp[which(y==0)] <- -mu[which(y==0)]
  }

  if (verbose == TRUE) {
    print("converged")
  }

  ## elements to return
  k   <- ncol(matrix(x))
  n   <- length(selectobs)
  reg$mu  <- mu
  reg$deviance <- -2 * sum(temp)/n
  reg$bic <- deviance + k*log(n)/n

  rownames(reg$coefficients) <- xnames

  # k = number of elements in x here
  # BIC would be BIC = deviance + k * ln(n)

  #returnlist <- list("coefficients" = b, "mu" = mu, "bic" = bic, "deviance" = deviance)
  if (saveX==TRUE) {
    reg[["x_resid"]] <- x_resid
    reg[["z_resid"]] <- z_resid
  }
  if (vcv) {
    if(!is.null(cluster)) {
      nclusters  <- nlevels(droplevels(cluster, exclude = if(anyNA(levels(cluster))) NULL else NA))
      het_matrix <- (1/nclusters)*cluster_matrix((y-mu)/sum(sqrt(mu)),cluster, x_resid)
      W          <- (1/nclusters)*(t(mu*x_resid) %*% x_resid)/sum(sqrt(mu))
      R <- try(chol(W), silent = FALSE)
      V          <- (1/nclusters)*chol2inv(R)%*%het_matrix%*%chol2inv(R)
      #V          <- (1/nclusters)*solve(W)%*%het_matrix%*%solve(W)
      V          <- nclusters / (nclusters-1) * V
    } else {
      e = y-mu
      het_matrix = (1/n) * t(x_resid*e)  %*% (x_resid*e)
      W          = (1/n) * (t(mu*x_resid) %*% x_resid)
      R          = try(chol(W), silent = TRUE)
      V          = (1/n) * chol2inv(R)%*%het_matrix%*%chol2inv(R)
      V          = (n/(n-1)) * V
    }
  }
  reg[["se"]] <- (diag(sqrt(V)))
  return(reg)
}
