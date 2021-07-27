#' General Penalized PPML Estimation
#'
#' \code{mlfitppml} is a general wrapper function for penalized PPML estimation, calling
#' \code{penhdfeppml}, \code{penhdfeppml_cluster} and \code{hdfeppml} as needed.
#'
#' @param lambdas Vector of penalty parameters.
#' @param IDs TODO: check what this does (probably for cross validation).
#' @param xval Logical. If \code{TRUE}, carries out cross-validation.
#' @param K TODO: check what this does.
#' @param vcv TODO: check this.
#' @inheritParams penhdfeppml
#'
#' @return A list (TODO: complete this).
#' @export
#'
#' @examples # TODO: add examples here.

mlfitppml = function(y,x,fes,lambdas,penalty=c("lasso","ridge"),tol=1e-8,hdfetol=1e-4,colcheck=TRUE,post=TRUE,cluster=NULL,
                     method=c("bic","iterative"),IDs=1:n,verbose=FALSE,
                     xval=FALSE,standardize=TRUE,vcv=TRUE,penweights=NULL,K=15) {

  set.seed(1)
  xnames <- colnames(x)
  n      <- length(y)

  # collinearity check: if option selected drop x's that are perfectly collinear beforehand
  if (colcheck==TRUE){
    include_x <- collinearity_check(y,x,fes,1e-6)
    x <- x[,include_x]
    colnames(x) <- xnames[include_x]
    xnames <- xnames[include_x]
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

  # method == "iterative" => use plugin method, which iterates on regressor-specific penalty weights.
  if (method == "iterative" & xval == FALSE) {

    # for storing current estimates
    pen_beta <- matrix(0,nrow = ncol(x), ncol = 1)

    # this is the subcommand that implements the iterative plugin estimator
    penreg   <- penhdfeppml_cluster(y=y,x=x,fes=fes,tol=tol,hdfetol=hdfetol,penalty=penalty,
                                    cluster=cluster,colcheck=FALSE,post=FALSE,verbose=verbose,K=K)

    ses <- matrix(NA,nrow = ncol(x), ncol = 1)

    # if "post", implement post-lasso PPML using selected covariates
    if (post) {
      x_select <- penreg$x_resid[,as.numeric(penreg$beta)!=0]

      # "pen_beta_pre" are raw lasso coefficients
      pen_beta_pre <- penreg$beta
      #x_select <- x[,as.numeric(penreg$beta)!=0]
      if(length(x_select)!=0){
        ppml_temp <- hdfeppml(y=y,x=x_select,fes=fes,tol=tol,hdfetol=hdfetol,mu=penreg$mu,
                              colcheck=FALSE,cluster=cluster)

        print(ppml_temp$coefficients)

        pen_beta[which(penreg$beta!=0),1]  <- ppml_temp$coefficients
        pen_bic   <- ppml_temp$bic
        print(ses[which(penreg$beta!=0),1])
        print(ppml_temp$se)
        ses[which(penreg$beta!=0),1] <- ppml_temp$se
      }
    }
    else{
      pen_beta <- penreg$beta
      pen_bic <- penreg$bic
    }
    pen_beta <- t(pen_beta)
    colnames(pen_beta) <- xnames

    # store results
    results <- list("beta" = t(pen_beta), "beta_pre" = t(pen_beta_pre), "deviance" = penreg$deviance, "bic" = penreg$bic, "lambda" = penreg$lambda, "phi" =penreg$phi, "ses" =t(ses))
  } else {

    # if method != "iterative", do the following:
    lambdas = sort(lambdas,decreasing=TRUE)
    pen_beta <- matrix(nrow = ncol(x), ncol = length(lambdas))
    pen_ses  <- matrix(nrow = ncol(x), ncol = length(lambdas))
    pen_bic  <- matrix(nrow = 1, ncol = length(lambdas))

    last_penbeta <- matrix(0,nrow=ncol(x),ncol=1)
    for (v in 1:length(lambdas)) {
      # compute lasso results for a single lambda

      print(lambdas[v])
      if (v==1) {
        penreg <- penhdfeppml(y=y,x=x,fes=fes,lambda=lambdas[v],tol=tol,hdfetol=hdfetol,
                              penalty=penalty,colcheck=FALSE,post=FALSE,standardize=standardize,method=method,cluster=cluster,penweights=penweights)
      } else {
        last_penbeta <- penreg$beta
        penreg <- penhdfeppml(y=y,x=penreg$x_resid,fes=fes,lambda=lambdas[v],tol=tol,hdfetol=hdfetol,
                              penalty=penalty,mu=penreg$mu,colcheck=FALSE,post=FALSE,standardize=standardize,method=method,cluster=cluster,penweights=penweights)
      }

      pen_beta[,v] <- penreg$beta

      # if "xval" is enabled, implement cross-validation (see xvalidate.R)
      if(xval==TRUE) {

        # opportunity to pass pen weights here.
        xval_reg   <- xvalidate(y=y,x=x,fes=fes,IDs=IDs,tol=tol,hdfetol=hdfetol,colcheck=TRUE,lambda=lambdas[v],cluster=cluster,
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
            pen_beta[,v] <- pen_beta[,v-1]
            pen_ses[,v]  <- pen_ses[,v-1]
            pen_bic[,v]  <- pen_bic[,v-1]
          } else{
            # pass mu, x, z as arguments here.
            if(length(x_select)!=0){
              ppml_temp <- hdfeppml(y = y, x = x_select, fes = fes, tol = tol, hdfetol = hdfetol,
                                    mu = penreg$mu, colcheck = FALSE, cluster = cluster, vcv = vcv)
              pen_beta_pre <- pen_beta

              pen_beta[which(penreg$beta!=0),v]  <- ppml_temp$coefficients
              pen_ses[which(penreg$beta!=0),v]   <- ppml_temp$se
              pen_bic[,v]   <- ppml_temp$bic
            }
          }
        }

      } else {
        pen_beta[,v] <- penreg$beta
        pen_bic[,v]  <- penreg$bic
      }
    }
    pen_beta <- t(pen_beta)
    colnames(pen_beta) <- xnames
    pen_beta_pre <- t(pen_beta_pre)
    colnames(pen_beta_pre) <- xnames
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

