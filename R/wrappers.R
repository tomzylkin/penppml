#' PPML Estimation for Data Frames
#'
#' These functions are thin wrappers around \code{mlfitppml}, \code{hdfeppml}, \code{penhdfeppml} and
#' \code{penhdfeppml_cluster}, providing a more convenient interface for data frames. Whereas the original
#' functions require some preliminary handling of data sets (\code{y} must be a vector, \code{x} must be a
#' matrix and \code{fes} must be provided in a list), the wrappers take a full data frame in the \code{data}
#' argument, and users can simply specify which variables correspond to y, x and the fixed effects, using
#' either variable names or column numbers.
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
#' @return The same as their respective parent functions.
#'
#' @examples
#' # To reduce run time, we keep only countries in the Americas:
#' americas <- countries$iso[countries$region == "Americas"]
#' # Now we can use our main functions on the reduced trade data set:
#' test1 <- mlfitppml2(data = trade[, -(5:6)],
#'                     dep = "export",
#'                     fixed = list(c("exp", "time"),
#'                                  c("imp", "time"),
#'                                  c("exp", "imp")),
#'                     selectobs = (trade$imp %in% americas) & (trade$exp %in% americas),
#'                     lambdas = c(0.01, 0.001, 0.0001))
#' test2 <- hdfeppml2(data = trade[, -(5:6)],
#'                    dep = "export",
#'                    fixed = list(c("exp", "time"),
#'                                 c("imp", "time"),
#'                                 c("exp", "imp")),
#'                    selectobs = (trade$imp %in% americas) & (trade$exp %in% americas))
#' test3 <- penhdfeppml2(data = trade[, -(5:6)],
#'                       dep = "export",
#'                       fixed = list(c("exp", "time"),
#'                                    c("imp", "time"),
#'                                    c("exp", "imp")),
#'                       lambda = 0.05,
#'                       selectobs = (trade$imp %in% americas) & (trade$exp %in% americas))
#' test4 <- penhdfeppml_cluster2(data = trade[, -(5:6)],
#'                               dep = "export",
#'                               fixed = list(c("exp", "time"),
#'                                            c("imp", "time"),
#'                                            c("exp", "imp")),
#'                               cluster = c("exp", "imp"),
#'                               selectobs = (trade$imp %in% americas) & (trade$exp %in% americas))
#'
#' @seealso [mlfitppml()]
#' @name wrappers

NULL

#' @export
#' @rdname wrappers
mlfitppml2 <- function(data, dep = 1, indep = NULL, fixed = NULL, cluster = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep,
                    fixed = fixed, cluster = cluster, selectobs = selectobs)
  # Final call to mlfitppml:
  mlfitppml(y = model$y, x = model$x, fes = model$fes, cluster = model$cluster, ...)
}

#' @export
#' @rdname wrappers
hdfeppml2 <- function(data, dep = 1, indep = NULL, fixed = NULL, cluster = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep,
                    fixed = fixed, cluster = cluster, selectobs = selectobs)
  # Final call to mlfitppml:
  hdfeppml(y = model$y, x = model$x, fes = model$fes, cluster = model$cluster, ...)
}

#' @export
#' @rdname wrappers
penhdfeppml2 <- function(data, dep = 1, indep = NULL, fixed = NULL, cluster = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep,
                    fixed = fixed, cluster = cluster, selectobs = selectobs)
  # Final call to mlfitppml:
  penhdfeppml(y = model$y, x = model$x, fes = model$fes, cluster = model$cluster, ...)
}

#' @export
#' @rdname wrappers
penhdfeppml_cluster2 <- function(data, dep = 1, indep = NULL, fixed = NULL, cluster = NULL, selectobs = NULL, ...) {
  # Initial call to genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep, fixed = fixed,
                    cluster = cluster, selectobs = selectobs)
  # Final call to mlfitppml:
  penhdfeppml_cluster(y = model$y, x = model$x, fes = model$fes, cluster = model$cluster, ...)
}
