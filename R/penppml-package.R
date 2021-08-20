#' penppml: Penalized Poisson Pseudo Maximum Likelihood Regression
#'
#' The `penppml` package is a set of tools that enables efficient estimation of penalized Poisson
#' Pseudo Maximum Likelihood (PPML) regressions, using lasso or ridge penalties, for models that
#' feature one or more sets of high-dimensional fixed effects (HDFE). The methodology is based on
#' Breinlich, Corradi, Rocha, Ruta, Santos Silva, and Zylkin (2021) and takes advantage of the method
#' of alternating projections of Gaure (2013) for dealing with HDFE, as well as the coordinate descent
#' algorithm of Friedman, Hastie and Tibshirani (2010) for fitting lasso regressions. The package is
#' also able to carry out cross-validation and to implement the plugin lasso of Belloni, Chernozhukov,
#' Hansen and Kozbur (2016).
#'
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
#' @docType package
#' @name penppml
NULL

## usethis namespace: start
#' @useDynLib penppml, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

