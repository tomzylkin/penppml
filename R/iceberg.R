#' Iceberg Lasso Implementation (in development)
#'
#' A function performs standard plugin lasso PPML estimation (without fixed effects) for several
#' dependent variables in a single step. This is still IN DEVELOPMENT: at the current stage, only
#' coefficient estimates are are provided and there is no support for clustered errors.
#'
#' This functions enables users to implement the "iceberg" step in the two-step procedure described in
#' Breinlich, Corradi, Rocha, Ruta, Santos Silva and Zylkin (2020). To do this after using the plugin
#' method in \code{mlfitppml}, just select all the variables with non-zero coefficients in
#' \code{dep} and the remaining regressors in \code{indep}. The function will then perform separate
#' lasso estimation on each of the selected dependent variables and report the coefficients.
#'
#' @param data A data frame containing all relevant variables.
#' @param dep A string with the names of the independent variables or their column numbers.
#' @param indep A vector with the names or column numbers of the regressors. If left unspecified,
#'              all remaining variables (excluding fixed effects) are included in the regressor matrix.
#' @param selectobs Optional. A vector indicating which observations to use (either a logical vector
#'     or a numeric vector with row numbers, as usual when subsetting in R).
#' @param ... Further arguments, including:
#' \itemize{
#' \item \code{tol}: Tolerance parameter for convergence of the IRLS algorithm.
#' \item \code{glmnettol}: Tolerance parameter to be passed on to \code{glmnet::glmnet}.
#' \item \code{penweights}: Optional: a vector of coefficient-specific penalties to use in plugin lasso.
#' \item \code{colcheck}: Logical. If \code{TRUE}, checks for perfect multicollinearity in \code{x}.
#' \item \code{K}: Maximum number of iterations.
#' \item \code{verbose}: Logical. If \code{TRUE}, prints information to the screen while evaluating.
#' \item \code{lambda}: Penalty parameter (a number).
#' \item \code{icepost}: Logical. If \code{TRUE}, it carries out a post-lasso estimation with just the
#'     selected variables and reports the coefficients from this regression.
#' }
#'
#' @return A matrix with coefficient estimates for all dependent variables.
#' @export
#'
#' @examples
#' \donttest{iceberg_results <- iceberg(data = trade[, -(1:6)],
#'                                     dep = c("ad_prov_14", "cp_prov_23", "tbt_prov_07",
#'                                             "tbt_prov_33", "tf_prov_41", "tf_prov_45"),
#'                                     selectobs = (trade$time == "2016"))}
#'
#' @inheritSection hdfeppml_int References

iceberg <- function(data, dep, indep = NULL, selectobs = NULL, ...) {
  # First we do the data handling with genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep, selectobs = selectobs)

  y_mat <- as.matrix(model$y)
  if(is.numeric(dep)){colnames(y_mat) <- colnames(data)[dep]} else {colnames(y_mat) <- dep}
  # Now we create the result matrix:
  iceberg_results <- matrix(NA, nrow = ncol(model$x), ncol = ncol(y_mat))
  rownames(iceberg_results) <- colnames(model$x)
  colnames(iceberg_results) <- colnames(y_mat)
  # Finally, we call plugin_lasso_int
  for (v in 1:ncol(y_mat)) {
    temp <- plugin_lasso_int(y = y_mat[, v], x = model$x, K = 15)
    iceberg_results[, v] <- temp$beta
  }
  return(iceberg_results)
}

#' Iceberg Lasso Implementation (in development)
#'
#' This is the internal function upon which the `iceberg` wrapper is built. It performs standard
#' plugin lasso PPML estimation without fixed effects, relying on \code{glmnet::glmnet}. As the other
#' internals in the package, it needs a y vector and an x matrix.
#'
#' @param y Dependent variable (a vector).
#' @param x Regressor matrix.
#' @param tol Tolerance parameter for convergence of the IRLS algorithm.
#' @param glmnettol Tolerance parameter to be passed on to \code{glmnet::glmnet}.
#' @param penweights Optional: a vector of coefficient-specific penalties to use in plugin lasso.
#' @param colcheck Logical. If \code{TRUE}, checks for perfect multicollinearity in \code{x}.
#' @param K Maximum number of iterations.
#' @param verbose Logical. If \code{TRUE}, prints information to the screen while evaluating.
#' @param lambda Penalty parameter (a number).
#' @param icepost Logical. If \code{TRUE}, it carries out a post-lasso estimation with just the
#'     selected variables and reports the coefficients from this regression.
#'
#' @return A list with 14 elements, including \code{beta}, which is the only one we use in the wrapper.
#' For a full list, see \link[glmnet]{glmnet}.

plugin_lasso_int <- function(y, x, tol = 1e-8,
                         glmnettol = 1e-12, penweights = NULL,
                         colcheck = FALSE, K = 50, verbose = FALSE, lambda = NULL, icepost = FALSE) {

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
#    print(iter)
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

    if (icepost == TRUE) {
      x_select <- x[, as.numeric(penreg$beta) != 0]
      if (length(x_select) != 0) {
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

  penreg[["beta"]] <- b
  penreg[["deviance"]] <- deviance
  penreg[["phi"]] <- phi
  penreg[["lambda"]] <- lambda / n

  return(penreg)
}
