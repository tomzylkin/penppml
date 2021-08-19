#' General Penalized PPML Estimation
#'
#' \code{mlfitppml} is a general-purpose wrapper function for penalized PPML estimation. This is a
#' flexible tool that allows users to select:
#' \itemize{
#'     \item Penalty type: either lasso or ridge.
#'     \item Penalty parameter: users can provide a single global value for lambda (a single regression
#'     is estimated), a vector of lambda values (the function estimates the regression using each of them,
#'     sequentially) or even coefficient-specific penalty weights.
#'     \item Method: plugin lasso estimates can be obtained directly from this function too.
#'     \item Cross-validation: if this option is enabled, the function uses IDs provided by the user
#'         to perform k-fold cross-validation and reports the resulting RMSE for all lambda values.
#'     }
#'
#' This function is a thin wrapper around \code{mlfitppml_int}, providing a more convenient interface for
#' data frames. Whereas the internal function requires some preliminary handling of data sets (\code{y}
#' must be a vector, \code{x} must be a matrix and \code{fes} must be provided in a list), the wrapper
#' takes a full data frame in the \code{data} argument, and users can simply specify which variables
#' correspond to y, x and the fixed effects, using either variable names or column numbers.
#'
#' For technical details on the algorithms used, see \link{hdfeppml} (post-lasso regression),
#' \link{penhdfeppml} (standard penalized regression), \link{penhdfeppml_cluster} (plugin lasso),
#' and \link{xvalidate} (cross-validation).
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
#' @param ... Further arguments, including:
#' \itemize{
#'     \item \code{penalty}: A string indicating the penalty type. Currently supported: "lasso" and "ridge".
#'     \item \code{method}: The user can set this equal to "plugin" to perform the plugin algorithm with
#'         coefficient-specific penalty weights (see details). Otherwise, a single global penalty is used.
#'     \item \code{post}: Logical. If \code{TRUE}, estimates a post-penalty regression with the
#'         selected variables.
#'     \item \code{xval}: Logical. If \code{TRUE}, cross-validation is performed using the IDs provided
#'         in the \code{IDs} argument as folds. Note that, by default, observations are assigned
#'         individual IDs, which makes the cross-validation algorithm very time-consuming.
#' }
#' For a full list of options, see \link{mlfitppml_int}.
#'
#' @inherit mlfitppml_int return
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
#'                     lambdas = c(0.01, 0.001),
#'                     tol = 1e-6, hdfetol = 1e-2)
#'
#' @inheritSection hdfeppml_int References

mlfitppml <- function(data, dep = 1, indep = NULL, fixed = NULL, cluster = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep,
                    fixed = fixed, cluster = cluster, selectobs = selectobs)
  # Final call to mlfitppml_int:
  mlfitppml_int(y = model$y, x = model$x, fes = model$fes, cluster = model$cluster, ...)
}


#' PPML Estimation with HDFE
#'
#' \code{hdfeppml} fits an (unpenalized) Poisson Pseudo Maximum Likelihood (PPML) model with
#' high-dimensional fixed effects (HDFE).
#'
#' This function is a thin wrapper around \link{hdfeppml_int}, providing a more convenient interface for
#' data frames. Whereas the internal function requires some preliminary handling of data sets (\code{y}
#' must be a vector, \code{x} must be a matrix and fixed effects \code{fes} must be provided in a list),
#' the wrapper takes a full data frame in the \code{data} argument, and users can simply specify which
#' variables correspond to y, x and the fixed effects, using either variable names or column numbers.
#'
#' More formally, \code{hdfeppml_int} performs iteratively re-weighted least squares (IRLS) on a
#' transformed model, as described in Correia, Guimarães and Zylkin (2020) and similar to the
#' \code{ppmlhdfe} package in Stata. In each iteration, the function calculates the transformed dependent
#' variable, partials out the fixed effects (calling \code{lfe::demeanlist}) and then solves a weighted
#' least squares problem (using fast C++ implementation).
#'
#' @inheritParams mlfitppml
#' @param ... Further options. For a full list, see \link{hdfeppml_int}.
#'
#' @inherit hdfeppml_int return
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
#'
#' @inheritSection hdfeppml_int References

hdfeppml <- function(data, dep = 1, indep = NULL, fixed = NULL, cluster = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep,
                    fixed = fixed, cluster = cluster, selectobs = selectobs)
  # Final call to mlfitppml_int:
  hdfeppml_int(y = model$y, x = model$x, fes = model$fes, cluster = model$cluster, ...)
}


#' One-Shot Penalized PPML Estimation with HDFE
#'
#' \code{penhdfeppml} fits a penalized PPML regression for a given type of penalty and a given
#' value of the penalty parameter.  The penalty can be either lasso or ridge, and the plugin method
#' can be enabled via the \code{method} argument.
#'
#' This function is a thin wrapper around \link{penhdfeppml_int}, providing a more convenient interface
#' for data frames. Whereas the internal function requires some preliminary handling of data sets (\code{y}
#' must be a vector, \code{x} must be a matrix and \code{fes} must be provided in a list), the wrapper
#' takes a full data frame in the \code{data} argument, and users can simply specify which variables
#' correspond to y, x and the fixed effects, using either variable names or column numbers.
#'
#' More formally, \code{penhdfeppml_int} performs iteratively re-weighted least squares (IRLS) on a
#' transformed model, as described in Breinlich, Corradi, Rocha, Ruta, Santos Silva and Zylkin (2021).
#' In each iteration, the function calculates the transformed dependent variable, partials out the fixed
#' effects (calling \code{lfe::demeanlist}) and then and then calls \code{glmnet::glmnet} if the selected
#' penalty is lasso (the default). If the user has selected ridge, the analytical solution is instead
#' computed directly using fast C++ implementation.
#'
#' For information on how the plugin lasso method works, see \link{penhdfeppml_cluster}.
#'
#' @inheritParams mlfitppml
#' @param ... Further options, including:
#' \itemize{
#'     \item \code{penalty}: A string indicating the penalty type. Currently supported: "lasso" and "ridge".
#'     \item \code{method}: The user can set this equal to "plugin" to perform the plugin algorithm with
#'         coefficient-specific penalty weights (see details). Otherwise, a single global penalty is used.
#' }
#' For a full list of options, see \link{penhdfeppml_int}.
#'
#' @inherit penhdfeppml_int return
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
#'
#' @inheritSection hdfeppml_int References

penhdfeppml <- function(data, dep = 1, indep = NULL, fixed = NULL, cluster = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep,
                    fixed = fixed, cluster = cluster, selectobs = selectobs)
  # Final call to mlfitppml_int:
  penhdfeppml_int(y = model$y, x = model$x, fes = model$fes, cluster = model$cluster, ...)
}


#' Plugin Lasso Estimation
#'
#' Performs plugin lasso - PPML estimation with HDFE. This is an internal function, called by
#' \code{mlfitppml} and \code{penhdfeppml} when users select the \code{method = "plugin"}
#' option, but it's made available as a stand-alone option for advanced users who may prefer to avoid
#' some overhead imposed by the wrappers.
#'
#' This function is a thin wrapper around \code{penppml_cluster_int}, providing a more convenient interface
#' for data frames. Whereas the internal function requires some preliminary handling of data sets (\code{y}
#' must be a vector, \code{x} must be a matrix and \code{fes} must be provided in a list), the wrapper
#' takes a full data frame in the \code{data} argument, and users can simply specify which variables
#' correspond to y, x and the fixed effects, using either variable names or column numbers.
#'
#' The plugin method uses coefficient-specific penalty weights that account for heteroskedasticity. The
#' penalty parameters are calculated automatically by the function using statistical theory - for a
#' brief discussion of this, see Breinlich, Corradi, Rocha, Ruta, Santos Silva and Zylkin (2021), and
#' for a more in-depth analysis, check Belloni, Chernozhukov, Hansen, and Kozbur (2016), which introduced
#' the specific implementation used in this package. Heuristically, the penalty parameters are set at
#' a level high enough so that the absolute value of the score for each regressor must be statistically
#' large relative to its standard error in order for the regressors to be selected.
#'
#' @inheritParams mlfitppml
#' @param cluster A string with the name of the clustering variable or a column number.
#'     It's also possible to input a vector with several variables, in which case the interaction of
#'     all of them is taken as the clustering variable. Note that this is NOT OPTIONAL in this case:
#'     our plugin algorithm requires clusters to be specified.
#' @param ... Further options. For a full list of options, see \link{penhdfeppml_cluster_int}.
#'
#' @inherit penhdfeppml_cluster_int return
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
#'
#' @inheritSection hdfeppml_int References

penhdfeppml_cluster <- function(data, dep = 1, indep = NULL, fixed = NULL, cluster = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep, fixed = fixed,
                    cluster = cluster, selectobs = selectobs)
  # Final call to mlfitppml_int:
  penhdfeppml_cluster_int(y = model$y, x = model$x, fes = model$fes, cluster = model$cluster, ...)
}
