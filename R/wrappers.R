#' General Penalized PPML Estimation
#'
#' \code{mlfitppml} is a general wrapper function for penalized PPML estimation.
#'
#' This function is a thin wrapper around \code{mlfitppml_int}, providing a more convenient interface for
#' data frames. Whereas the internal function requires some preliminary handling of data sets (\code{y}
#' must be a vector, \code{x} must be a matrix and \code{fes} must be provided in a list), the wrapper
#' takes a full data frame in the \code{data} argument, and users can simply specify which variables
#' correspond to y, x and the fixed effects, using either variable names or column numbers.
#'
#' @param data A data frame containing all relevant variables.
#' @param dep A string with the name of the independent variable or a column number.
#' @param indep A vector with the names or column numbers of the regressors. If left unspecified,
#'              all remaining variables (excluding fixed effects) are included in the regressor matrix.
#' @param fixed A vector with the names or column numbers of factor variables identifying the fixed effects,
#'     or a list with the desired interactions between variables in \code{data}.
#' @param selectobs Optional. A vector indicating which observations to use (either a logical vector
#'     or a numeric vector with row numbers, as usual when subsetting in R).
#' @param cluster Optional. A string with the name of the clustering variable or a column number.
#'     It's also possible to input a vector with several variables, in which case the interaction of
#'     all of them is taken as the clustering variable.
#' @param ... Further arguments, to be passed on to the main function.
#'
#' @return TODO: add this.
#' @export
#'
#' @examples
#' # To reduce run time, we keep only countries in the Americas:
#' americas <- countries$iso[countries$region == "Americas"]
#' # Now we can use our main functions on the reduced trade data set:
#' test <- mlfitppml(data = trade[, -(5:6)],
#'                     dep = "export",
#'                     fixed = list(c("exp", "time"),
#'                                  c("imp", "time"),
#'                                  c("exp", "imp")),
#'                     selectobs = (trade$imp %in% americas) & (trade$exp %in% americas),
#'                     lambdas = c(0.01, 0.001, 0.0001),
#'                     tol = 1e-5, hdfetol = 1e-1)

mlfitppml <- function(data, dep = 1, indep = NULL, fixed = NULL, cluster = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep,
                    fixed = fixed, cluster = cluster, selectobs = selectobs)
  # Final call to mlfitppml_int:
  mlfitppml_int(y = model$y, x = model$x, fes = model$fes, cluster = model$cluster, ...)
}


#' Poisson Pseudo Maximum Likelihood Estimation
#'
#' \code{hdfeppml} implements (unpenalized) PPML estimation in the presence of high-dimensional
#' fixed effects.
#'
#' This function is a thin wrapper around \code{hdfeppml_int}, providing a more convenient interface for
#' data frames. Whereas the internal function requires some preliminary handling of data sets (\code{y}
#' must be a vector, \code{x} must be a matrix and \code{fes} must be provided in a list), the wrapper
#' takes a full data frame in the \code{data} argument, and users can simply specify which variables
#' correspond to y, x and the fixed effects, using either variable names or column numbers.
#'
#' Internally, \code{hdfeppml_int} performs iteratively re-weighted least squares (IRLS) on a transformed
#' model, as described in Breinlich, Corradi, Rocha, Ruta, Santos Silva and Zylkin (2021). In each
#' iteration, the function calculates the transformed dependent variable, partials out the fixed effects
#' (calling \code{lfe::demeanlist}) and then solves a weighted least squares problem (using fast C++
#' implementation).
#'
#' @inheritParams mlfitppml
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{coefficients}: coefficient (beta) estimates.
#'   \item \code{residuals}: residuals of the model.
#'   \item \code{mu}: (TODO: check this)
#'   \item \code{deviance}: (TODO: check this)
#'   \item \code{bic}: (TODO: check this)
#'   \item \code{x_resid}: matrix of demeaned regressors.
#'   \item \code{z_resid}: vector of demeaned (transformed) dependent variable.
#'   \item \code{se}: standard errors of the coefficients.
#' }
#' @export
#'
#' @examples
#' # To reduce run time, we keep only countries in the Americas:
#' americas <- countries$iso[countries$region == "Americas"]
#' test <- hdfeppml(data = trade[, -(5:6)],
#'                    dep = "export",
#'                    fixed = list(c("exp", "time"),
#'                                 c("imp", "time"),
#'                                 c("exp", "imp")),
#'                    selectobs = (trade$imp %in% americas) & (trade$exp %in% americas))

hdfeppml <- function(data, dep = 1, indep = NULL, fixed = NULL, cluster = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep,
                    fixed = fixed, cluster = cluster, selectobs = selectobs)
  # Final call to mlfitppml_int:
  hdfeppml_int(y = model$y, x = model$x, fes = model$fes, cluster = model$cluster, ...)
}


#' One-Shot Penalized PPML Estimation
#'
#' \code{penhdfeppml} computes a penalized PPML model for a given type of penalty and a given
#' value of the penalty parameter.
#'
#' This function is a thin wrapper around \code{penppml_int}, providing a more convenient interface
#' for data frames. Whereas the internal function requires some preliminary handling of data sets (\code{y}
#' must be a vector, \code{x} must be a matrix and \code{fes} must be provided in a list), the wrapper
#' takes a full data frame in the \code{data} argument, and users can simply specify which variables
#' correspond to y, x and the fixed effects, using either variable names or column numbers.
#'
#' @inheritParams mlfitppml
#'
#' @return TODO: add this.
#' @export
#'
#' @examples
#' # To reduce run time, we keep only countries in the Americas:
#' americas <- countries$iso[countries$region == "Americas"]
#' test <- penhdfeppml(data = trade[, -(5:6)],
#'                       dep = "export",
#'                       fixed = list(c("exp", "time"),
#'                                    c("imp", "time"),
#'                                    c("exp", "imp")),
#'                       lambda = 0.05,
#'                       selectobs = (trade$imp %in% americas) & (trade$exp %in% americas))

penhdfeppml <- function(data, dep = 1, indep = NULL, fixed = NULL, cluster = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep,
                    fixed = fixed, cluster = cluster, selectobs = selectobs)
  # Final call to mlfitppml_int:
  penhdfeppml_int(y = model$y, x = model$x, fes = model$fes, cluster = model$cluster, ...)
}


#' Plugin Lasso Estimation
#'
#' Estimates coefficient-specific penalty weights that account for heteroskedasticity.
#' Called by \code{mlfitppml} and \code{penhdfeppml}.
#'
#' This function is a thin wrapper around \code{penppml_cluster_int}, providing a more convenient interface
#' for data frames. Whereas the internal function requires some preliminary handling of data sets (\code{y}
#' must be a vector, \code{x} must be a matrix and \code{fes} must be provided in a list), the wrapper
#' takes a full data frame in the \code{data} argument, and users can simply specify which variables
#' correspond to y, x and the fixed effects, using either variable names or column numbers.
#'
#' @inheritParams mlfitppml
#'
#' @return TODO: add this.
#' @export
#'
#' @examples
#' # To reduce run time, we keep only countries in the Americas:
#' americas <- countries$iso[countries$region == "Americas"]
#' test <- penhdfeppml_cluster(data = trade[, -(5:6)],
#'                               dep = "export",
#'                               fixed = list(c("exp", "time"),
#'                                            c("imp", "time"),
#'                                            c("exp", "imp")),
#'                               cluster = c("exp", "imp"),
#'                               selectobs = (trade$imp %in% americas) & (trade$exp %in% americas),
#'                               tol = 1e-5, hdfetol = 1e-1)

penhdfeppml_cluster <- function(data, dep = 1, indep = NULL, fixed = NULL, cluster = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep, fixed = fixed,
                    cluster = cluster, selectobs = selectobs)
  # Final call to mlfitppml_int:
  penhdfeppml_cluster_int(y = model$y, x = model$x, fes = model$fes, cluster = model$cluster, ...)
}
