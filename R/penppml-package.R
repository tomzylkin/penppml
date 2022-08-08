#' @section Functions:
#' The workhorse of this package is the \code{mlfitppml} function, which allows users to carry out
#' penalized HDFE-PPML estimation with a wide variety of options. The syntax is very simple, allowing
#' users to select a data frame with all the relevant variables and then select dependent, independent
#' and fixed effects variables by name or column number.
#'
#' In addition, the internals \code{hdfeppml} (post-lasso regression), \code{penhdfeppml} (penalized
#' regression for a single lambda), \code{penhdfeppml_cluster} (plugin lasso), and \code{xvalidate} (cross-
#' validation) are made available on a stand-alone basis for advanced users.
#'
#' The package also includes alternative versions of \code{mlfitppml}, \code{hdfeppml}, \code{penhdfeppml}
#' and \code{penhdfeppml_cluster}. These (\code{mlfitppml_int}, \code{hdfeppml_int}, \code{penhdfeppml_int}
#' and \code{penhdfeppml_cluster_int}) use an alternative syntax: users must provide the dependent variable
#' in a vector, the regressors in a matrix and the fixed effects in a list.
#'
#' Finally, support for the iceberg lasso method in Breinlich, Corradi, Rocha, Ruta, Santos Silva,
#' and Zylkin (2021) is in development and can be accessed at its current stage via the \code{iceberg}
#' function.
#'
#' @inheritSection hdfeppml_int References
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

## usethis namespace: start
#' @useDynLib penppml, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom glmnet glmnet
## usethis namespace: end
NULL
