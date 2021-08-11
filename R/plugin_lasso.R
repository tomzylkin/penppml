#' Iceberg Lasso Implementation (in development)
#'
#' TODO.
#'
#' @param data A data frame containing all relevant variables.
#' @param dep A string with the name of the independent variable or a column number.
#' @param indep A vector with the names or column numbers of the regressors. If left unspecified,
#'              all remaining variables (excluding fixed effects) are included in the regressor matrix.
#' @param selectobs Optional. A vector indicating which observations to use (either a logical vector
#'     or a numeric vector with row numbers, as usual when subsetting in R).
#' @param ... Further arguments, to be passed on to the main function.
#'
#'
#' @return A list with (TODO).
#' @export
#'
#' @examples # TODO

plugin_lasso <- function(data, dep = 1, indep = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep, selectobs = selectobs)
  # Final call to plugin_lasso_int:
  penhdfeppml_cluster_int(y = model$y, x = model$x, ...)
}



#' Iceberg Lasso Implementation (in development)
#'
#' TODO.
#'
#' @param y Dependent variable (a vector).
#' @param x Regressor matrix.
#' @param tol Tolerance parameter.
#' @param glmnettol Tolerance parameter to be passed on to glmnet.
#' @param penweights TODO: check what this does (parameter-specific penalties in plugin lasso?).
#' @param colcheck Logical. If \code{TRUE}, checks for perfect multicollinearity in \code{x}.
#' @param K TODO: check what this does (number of iterations for plugin lasso, probably).
#' @param verbose If TRUE, prints information to the screen while evaluating.
#' @param lambda Penalty parameter (a number).
#' @param phipost TODO: check this (something with glmnet).
#'
#' @return A list with (TODO).
#' @export
#'
#' @examples # TODO.

plugin_lasso_int <- function(y, x, tol = 1e-8,
                         glmnettol = 1e-12, penweights = NULL,
                         colcheck = FALSE, K = 50, verbose = FALSE, lambda = NULL, phipost = FALSE) {

  x <- data.matrix(x)
  y <- as.matrix(y)
  xnames <- colnames(x)
  n <- length(y)
  k <- ncol(x)

  if (is.null(lambda)) {
    c <- 1.1
    gamma <- 0.1 / log(n)
    lambda <- c * sqrt(n) * stats::qnorm(1 - gamma / (2 * k))
  }

  b <- matrix(NA, nrow = ncol(x), ncol = 1)  # fix this later.
  rownames(b) <- colnames(x)
  include_x <- 1:ncol(x)

  if (colcheck == TRUE){
    w <- matrix(1/n, n)
    w <- w[1:length(w)]
    check <- stats::lm.wfit(x, y, w) #faster than lmfit
    check$coefficients
    include_x <- which(!is.na(check$coefficients))
    x <- x[, include_x]
  }

  # number of obs (needed for deviance)
  n <- length(y)

  ## estimation algorithm
  crit <- 1
  old_deviance <- 0
  iter <- 0

  while (crit > tol & iter < K) {
    iter <- iter + 1
    print(iter)
    if (iter == 1) {
      e <- y - mean(y)
    }
    e <- as.numeric(e)

    het_matrix <- (1/n) * t(x * e)  %*% (x * e)

    phi <- sqrt(diag(het_matrix))

    lambda_glmnet <- lambda / n * sum(phi) / k

    if (verbose == TRUE) {
      print(phi)
    }

    penreg <- glmnet::glmnet(x = x, y = y, weights = rep(1/n, n), lambda = lambda_glmnet, thresh = glmnettol,
                     penalty.factor = phi, standardize = FALSE)

    if (phipost == TRUE){
      x_select <- x[, as.numeric(penreg$beta) != 0]
      if(length(x_select) != 0){
        b_temp <- rep(0, length(include_x))
        b_temp[as.numeric(penreg$beta) != 0] <- fastolsCpp(x_select, y)
        b[include_x] <- b_temp
      }
    }
    else{
      b[include_x] <- penreg$beta
    }
    e <- (y - mean(y)) - (x - colMeans(x)) %*% b[include_x]

    # calculate deviance
    temp <-  -(e^2)

    deviance <- -2 * sum(temp)/n

    if (deviance < 0) deviance <- 0

    delta_deviance <- old_deviance - deviance

    if (!is.na(delta_deviance) & (deviance < 0.1 * delta_deviance)) {
      delta_deviance = deviance
    }

    denom_crit <- max(c(min(c(deviance, old_deviance)), 0.1))
    crit <- abs(delta_deviance) / denom_crit

    #print(deviance)

    old_deviance <- deviance
  }
  ## elements to return
  k   <- ncol(matrix(x))
  n   <- length(y)
  select_x <- which(b != 0)

  print(k)
  print(b)

  penreg[["beta"]] <- b
  penreg[["deviance"]] <- deviance
  penreg[["phi"]] <- phi
  penreg[["lambda"]] <- lambda / n

  return(penreg)
}
