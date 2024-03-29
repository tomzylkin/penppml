% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fastridgeCpp}
\alias{fastridgeCpp}
\title{Faster Ridge Regression}
\usage{
fastridgeCpp(X, y, lambda)
}
\arguments{
\item{X}{Regressor matrix.}

\item{y}{Dependent variable (a vector).}

\item{lambda}{Penalty parameter (a number).}
}
\value{
The vector of parameter (beta) estimates.
}
\description{
Finds Ridge solutions using C++.
}
