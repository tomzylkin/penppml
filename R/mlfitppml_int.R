#' General Penalized PPML Estimation
#'
#' \code{mlfitppml_int} is the internal wrapper called by \code{mlfitppml} for penalized PPML estimation.
#' This in turn calls \code{penhdfeppml_int}, \code{penhdfeppml_cluster_int} and \code{hdfeppml_int}
#' as needed. It takes a vector with the dependent variable, a regressor matrix and a set of fixed
#' effects (in list form: each element in the list should be a separate HDFE). This is a flexible tool
#' that allows users to select:
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
#' For technical details on the algorithms used, see \link{hdfeppml_int} (post-lasso regression),
#' \link{penhdfeppml_int} (standard penalized regression), \link{penhdfeppml_cluster_int} (plugin lasso),
#' and \link{xvalidate} (cross-validation).
#'
#' @param lambdas Vector of penalty parameters.
#' @param IDs A vector of fold IDs for k-fold cross validation. If left unspecified, each observation
#'    is assigned to a different fold (warning: this is likely to be very resource-intensive).
#' @param xval Logical. If \code{TRUE}, it carries out cross-validation.
#' @param K Maximum number of iterations for the plugin algorithm to converge.
#' @param vcv Logical. If \code{TRUE} (the default), the post-estimation model includes standard errors.
#' @param phipost Logical. If \code{TRUE}, the plugin coefficient-specific penalty weights are iteratively
#' calculated using estimates from a post-penalty regression when \code{method == "plugin"}. Otherwise,
#' these are calculated using estimates from a penalty regression.
#' @param colcheck_x Logical. If \code{TRUE}, this checks collinearity between the independent variables and drops the
#' collinear variables.
#' @param colcheck_x_fes Logical. If \code{TRUE}, this checks whether the independent variables are perfectly explained
#' by the fixed effects drops those that are perfectly explained.
#' @inheritParams penhdfeppml_int
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{beta}: if \code{post = FALSE}, a \code{length(lambdas)} x \code{ncol(x)} matrix with
#'       coefficient (beta) estimates from the penalized regressions. If \code{post = TRUE}, this is
#'       the matrix of coefficients from the post-penalty regressions.
#'   \item \code{beta_pre}: if \code{post = TRUE}, a \code{length(lambdas)} x \code{ncol(x)} matrix with
#'       coefficient (beta) estimates from the penalized regressions.
#'   \item \code{bic}: Bayesian Information Criterion.
#'   \item \code{lambdas}: vector of penalty parameters.
#'   \item \code{ses}: standard errors of the coefficients of the post-penalty regression. Note that
#'       these are only provided when \code{post = TRUE}.
#'   \item \code{rmse}: if \code{xval = TRUE}, a matrix with the root mean squared error (RMSE - column 2)
#'       for each value of lambda (column 1), obtained by cross-validation.
#'   \item \code{phi}: coefficient-specific penalty weights (only if \code{method == "plugin"}).
#' }
#' @export
#'
#' @examples
#' # First, we need to transform the data (this is what mlfitppml handles internally). Start by
#' # filtering the data set to keep only countries in the Americas:
#' americas <- countries$iso[countries$region == "Americas"]
#' trade <- trade[(trade$imp %in% americas) & (trade$exp %in% americas), ]
#' # Now generate the needed x, y and fes objects:
#' y <- trade$export
#' x <- data.matrix(trade[, -1:-6])
#' fes <- list(exp_time = interaction(trade$exp, trade$time),
#'             imp_time = interaction(trade$imp, trade$time),
#'             pair     = interaction(trade$exp, trade$imp))
#' # Finally, we try mlfitppml_int with a lasso penalty (the default) and two lambda values:
#' reg <- mlfitppml_int(y = y, x = x, fes = fes, lambdas = c(0.1, 0.01))
#'
#' # We can also try plugin lasso:
#' \donttest{reg <- mlfitppml_int(y = y, x = x, fes = fes, cluster = fes$pair, method = "plugin")}
#'
#' # For an example with cross-validation, please see the vignette.
#'
#' @inheritSection hdfeppml_int References

mlfitppml_int = function(y, x, fes, lambdas, penalty = "lasso", tol = 1e-8, hdfetol = 1e-4, colcheck_x = FALSE, colcheck_x_fes = TRUE,
                     post = TRUE, cluster = NULL, method = "bic", IDs = 1:n, verbose = FALSE, xval = FALSE,
                     standardize = TRUE, vcv = TRUE, phipost=TRUE, penweights = NULL, K = 15, gamma_val=NULL, mu=NULL) {

  xnames <- colnames(x)
  n      <- length(y)

  # collinearity check: if option selected drop x's that are perfectly collinear beforehand
  if (colcheck_x==TRUE & colcheck_x_fes==TRUE){
    include_x <- collinearity_check(y,x,fes,1e-6, colcheck_x=colcheck_x, colcheck_x_fes=colcheck_x_fes)
    x <- x[,include_x]
    colnames(x) <- xnames[include_x]
    xnames <- xnames[include_x]
    colcheck_x_post = FALSE
    colcheck_x_fes_post = FALSE
  }
  if (colcheck_x==FALSE & colcheck_x_fes==FALSE){
    colcheck_x_post = TRUE
    colcheck_x_fes_post = TRUE
  }
  if (colcheck_x==TRUE & colcheck_x_fes==FALSE){
    include_x <- collinearity_check(y,x,fes,1e-6, colcheck_x=colcheck_x, colcheck_x_fes=colcheck_x_fes)
    x <- x[,include_x]
    colnames(x) <- xnames[include_x]
    xnames <- xnames[include_x]
    colcheck_x_post = TRUE
    colcheck_x_fes_post = FALSE
  }
  if (colcheck_x==FALSE & colcheck_x_fes==TRUE){
    include_x <- collinearity_check(y,x,fes,1e-6, colcheck_x=colcheck_x, colcheck_x_fes=colcheck_x_fes)
    x <- x[,include_x]
    colnames(x) <- xnames[include_x]
    xnames <- xnames[include_x]
    colcheck_x_post = FALSE
    colcheck_x_fes_post = TRUE
  }


  # "xval" = "cross-validation". If this option is selected, set up vectors to store cross-validation results
  if (xval==TRUE) {
    xval_rmse  <- matrix(nrow = 1, ncol = length(lambdas))
    xval_dev   <- matrix(nrow = 1, ncol = length(lambdas))
  }
  # if "ridge" regression is selected there is no need to post estimates
  if (penalty=="ridge") {
    post=FALSE
  }

  # method == "plugin" => use plugin method, which iterates on regressor-specific penalty weights.
  if (method == "plugin" & xval == FALSE) {
    if(is.null(gamma_val)){gamma_val <- 0.1/log(n)}
    # for storing current estimates
    pen_beta <- matrix(0,nrow = ncol(x), ncol = 1)

    # this is the subcommand that implements the iterative plugin estimator
    penreg   <- penhdfeppml_cluster_int(y=y,x=x,fes=fes,tol=tol,hdfetol=hdfetol,penalty=penalty,
                                    cluster=cluster,colcheck_x=FALSE,colcheck_x_fes=FALSE,post=FALSE,verbose=verbose,phipost=phipost,K=K, gamma_val=gamma_val, mu=mu)

    ses <- matrix(NA,nrow = ncol(x), ncol = 1)

    # if "post", implement post-lasso PPML using selected covariates
    if (post) {
      x_select <- penreg$x_resid[,as.numeric(penreg$beta)!=0]
      # "pen_beta_pre" are raw lasso coefficients
      pen_beta_pre <- penreg$beta
      #x_select <- x[,as.numeric(penreg$beta)!=0]
      if(length(x_select)!=0){
        ppml_temp <- hdfeppml_int(y=y,x=x_select,fes=fes,tol=tol,hdfetol=hdfetol,mu=penreg$mu,
                              colcheck_x=colcheck_x_post, colcheck_x_fes = colcheck_x_fes_post, cluster=cluster)

        pen_beta[which(penreg$beta!=0),1]  <- ppml_temp$coefficients
        pen_bic   <- ppml_temp$bic
        if(verbose==TRUE){print(ppml_temp$se)}
        ses[which(penreg$beta!=0),1] <- ppml_temp$se
      }
    } else {
      pen_beta <- penreg$beta
      pen_bic <- penreg$bic
    }
    pen_beta <- t(pen_beta)
    colnames(pen_beta) <- xnames
    # The following reproduces lines 208-211 for plugin lasso.
    if (post) {
      pen_beta_pre <- t(pen_beta_pre)
      colnames(pen_beta_pre) <- xnames
    }
    # store results, conditionally, because post may not be true
    if(post){results <- list("beta" = t(pen_beta), "beta_pre" = t(pen_beta_pre), "deviance" = penreg$deviance, "bic" = penreg$bic, "lambda" = penreg$lambda, "phi" =penreg$phi, "ses" =t(ses))} else {
      results <- list("beta" = t(pen_beta), "deviance" = penreg$deviance, "bic" = penreg$bic, "lambda" = penreg$lambda, "phi" =penreg$phi, "ses" =t(ses))
    }
  } else {

    # if method != "plugin", do the following:
    lambdas = sort(lambdas,decreasing=TRUE)
    pen_beta <- matrix(nrow = ncol(x), ncol = length(lambdas))
    if (post==TRUE) pen_beta_pre <- matrix(nrow = ncol(x), ncol = length(lambdas))
    pen_ses  <- matrix(nrow = ncol(x), ncol = length(lambdas))
    pen_bic  <- matrix(nrow = 1, ncol = length(lambdas))

    last_penbeta <- matrix(0,nrow=ncol(x),ncol=1)
    for (v in 1:length(lambdas)) {
      # compute lasso results for a single lambda
      print(lambdas[v])
      if (v==1) {
        penreg <- penhdfeppml_int(y=y,x=x,fes=fes,lambda=lambdas[v],tol=tol,hdfetol=hdfetol,
                              penalty=penalty,colcheck_x=FALSE,colcheck_x_fes=FALSE,post=FALSE,standardize=standardize,method=method,cluster=cluster,penweights=penweights, mu=mu)
      } else {
        last_penbeta <- penreg$beta
        penreg <- penhdfeppml_int(y=y,x=penreg$x_resid,fes=fes,lambda=lambdas[v],tol=tol,hdfetol=hdfetol,
                              penalty=penalty,mu=penreg$mu,colcheck_x=FALSE,colcheck_x_fes=FALSE,post=FALSE,standardize=standardize,method=method,cluster=cluster,penweights=penweights)
      }
      pen_beta[,v] <- penreg$beta

      # if "xval" is enabled, implement cross-validation (see xvalidate.R)
      if(xval==TRUE) {
        # opportunity to pass pen weights here.
        xval_reg   <- xvalidate(y=y,x=x,fes=fes,IDs=IDs,tol=tol,hdfetol=hdfetol,colcheck_x=FALSE,colcheck_x_fes=FALSE,lambda=lambdas[v],cluster=cluster,
                                init_mu=penreg$mu,init_x=x,init_z=penreg$z_resid,verbose=verbose,standardize=standardize,penalty=penalty,method=method,penweights=penweights)

        xval_rmse[,v] <- xval_reg$rmse
      }

      x_select <- penreg$x_resid[,as.numeric(penreg$beta)!=0]

      ##x_select <- x[,penreg$beta!=0]

      # compute post-lasso estimates if option enabled
      if (post) {
        if (sum(penreg$beta!=0)>0) {

          same <- min(((penreg$beta==0)==(last_penbeta==0)))

          # if selected x's same as for previous value of lambda, carry forward result.
          if (same==1 & v>1) {
            pen_beta_pre[,v] <- penreg$beta
            pen_beta[,v] <- pen_beta[,v-1]
            pen_ses[,v]  <- pen_ses[,v-1]
            pen_bic[,v]  <- pen_bic[,v-1]
          } else {
            # pass mu, x, z as arguments here.
            if(length(x_select)!=0){
              ppml_temp <- hdfeppml_int(y = y, x = x_select, fes = fes, tol = tol, hdfetol = hdfetol,
                                    mu = penreg$mu, colcheck_x_fes = colcheck_x_fes_post, colcheck_x = colcheck_x_post, cluster = cluster, vcv = vcv)
              pen_beta_pre[,v] <- penreg$beta

              pen_beta[which(penreg$beta!=0),v]  <- ppml_temp$coefficients
              pen_ses[which(penreg$beta!=0),v]   <- ppml_temp$se
              pen_bic[,v]   <- ppml_temp$bic
            }
          }
        } else pen_beta_pre[,v] <- penreg$beta
          # I've added the following to solve a bug: when post is TRUE and no variables are selected,
          # the beta_pre object should still be overwritten by 0s. Otherwise, NAs remain.
      } else {
        pen_beta[,v] <- penreg$beta
        pen_bic[,v]  <- penreg$bic
      }
    }
    pen_beta <- t(pen_beta)
    colnames(pen_beta) <- xnames
    if (post) {
      pen_beta_pre <- t(pen_beta_pre)
      colnames(pen_beta_pre) <- xnames
    }
    bic <- cbind(lambdas,t(pen_bic))
    colnames(bic) <- cbind("lambda","bic")


    # return results
    if (post) {
      #colnames(pen_beta_pre) <- xnames
      results <- list("beta" = t(pen_beta), "beta_pre" = t(pen_beta_pre), "bic" =  bic, "lambdas" = lambdas, "ses" = pen_ses)
    }
    else {
      results <- list("beta" = t(pen_beta), "bic" =  bic, "lambdas" = lambdas, "ses" = pen_ses)
    }
    if (xval) {
      results[["rmse"]]     <- cbind(lambdas,t(xval_rmse))
      #results[["deviance"]] <- cbind(lambdas,t(xval_dev))
    }
  }
  return(results)
}

