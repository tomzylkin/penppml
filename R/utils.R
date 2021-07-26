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


#' Checking for Perfect Multicollinearity
#'
#' \code{collinearity_check} checks for perfect multicollinearity in a model with high-dimensional
#' fixed effects.
#'
#' @param y Numeric vector with the dependent variable.
#' @param x Matrix whose columns will be checked for multicollinearity.
#' @param fes List of factors defining the fixed effects in the model.
#' @param hdfetol A number, passed on to \code{lfe::demeanlist} as a tolerance for the centering.
#' @param selectobs Currently unused.
#'
#' @return A numeric vector containing the variables that pass the collinearity check.
#' @export
#'
#' @examples
#' y <- trade$export
#' x <- data.matrix(trade[,-1:-9])
#' fes <- genfes(trade,
#'               f1 = c("exp",  "imp", "exp"),
#'               f2 = c("time",  "time", "imp"))
#'
#' collinearity_check(y, x, fes, hdfetol = 1e-6)

collinearity_check <- function(y, x, fes, hdfetol, selectobs = NULL) {
  mu  <- (y + mean(y)) / 2
  z   <- (y - mu) / mu + log(mu)
  reg_z  <- matrix(z)
  reg_x  <- x
  mu  <- (y + mean(y)) / 2

  z_resid <- lfe::demeanlist(reg_z, fes, weights = sqrt(mu), eps = hdfetol)
  x_resid <- lfe::demeanlist(reg_x, fes, weights = sqrt(mu), eps = hdfetol)

  check <- stats::lm.wfit(x_resid, z_resid, mu)
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

cluster_matrix <- function(e, cluster, x) {
  K <- ncol(x)
  vars      <- data.frame(e = e, cluster = factor(cluster, exclude = TRUE), x = x)
  vars      <- vars[order(vars$cluster),]
  vars$indx <- with(vars, ave(seq_along(cluster), cluster, FUN = seq_along))
  vars      <- tidyr::complete(vars, cluster, indx, fill = list(e = 0, x = 0))
  vars      <- vars %>% tidyr::fill(tidyr::everything(), 0)

  T <- max(vars$indx)

  if(is.null(K)){
    X <- data.matrix(vars[, 4])
  } else{
    X <- data.matrix(vars[, 4:(3 + K)])
  }

  e <- as.matrix(vars[, 3])
  rm(vars)

  ee   <- manyouter(e, e, T)  #many outer products for vectors of length T
  XeeX <- xeex(X, ee)        #XeeX matrix for cluster-robust SEs
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


#' Finding Ridge Regression Solutions
#'
#' A wrapper around \code{fastridgeCpp}, for faster computation of the analytical solution for
#' ridge regression.
#'
#' @param x Regressor matrix.
#' @param y Dependent variable (a numeric vector).
#' @param weights Vector of weights.
#' @param lambda Penalty parameter.
#' @param standardize Logical. If \code{TRUE}, x is standardized using the \code{weights}.
#'
#' @return
#' @export
#'
#' @examples

fastridge <- function(x, y, weights = rep(1/n, n), lambda, standardize = TRUE){
  n <- length(y)
  if (standardize) {
    b <- fastridgeCpp(sqrt(weights) * standardize_wt(x, weights), sqrt(weights) * y, lambda)
    beta <- b * faststddev(x, weights)
  } else {
    beta <- fastridgeCpp(sqrt(weights) * x, sqrt(weights) * y, lambda)
  }
  result <- list("beta" = beta)
  return(result)
}


#' Generating a List of Fixed Effects
#'
#' \code{genfes} generates a list of fixed effects by creating interactions of paired factors.
#'
#' @param data A data frame including the factors.
#' @param f1 A character vector with the names of variables to serve as first term in the interaction.
#' @param f2 A character vector with the names of variables (second term).
#'
#' @return A list containing the interacted factors.
#' @export
#'
#' @examples
#' data <- data.frame(exporter = c("Spain", "Spain", "Portugal", "France"),
#'                    importer = c("China", "Swaziland", "Mozambique", "Tunisia"),
#'                    year = c("1998", "1999", "1998", "1999"))
#' genfes(data, f1 = c("exporter",  "importer", "exporter"), f2 = c("year",  "year", "importer"))

genfes <- function(data, f1, f2) {
  fes <- list()

  if (length(f1) != length(f2)) {
    cat("ERROR: the length of f1 and f2 should be the same.")
    stop()
  }

  for (i in 1:length(f1)) {
    fes[[paste(f1[i], f2[i], sep = "_")]] <- interaction(data[[f1[i]]], data[[f2[i]]])
  }
  return(fes)
}

