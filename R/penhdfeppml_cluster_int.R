#' Plugin Lasso Estimation
#'
#' Performs plugin lasso - PPML estimation with HDFE. This is an internal function, called by \code{mlfitppml_int} and
#' \code{penhdfeppml_int} when users select the \code{method = "plugin"} option, but it's made available
#' as a stand-alone option for advanced users who may prefer to avoid some overhead imposed by the
#' wrappers.
#'
#' The plugin method uses coefficient-specific penalty weights that account for heteroskedasticity. The
#' penalty parameters are calculated automatically by the function using statistical theory - for a
#' brief discussion of this, see Breinlich, Corradi, Rocha, Ruta, Santos Silva and Zylkin (2021), and
#' for a more in-depth analysis, check Belloni, Chernozhukov, Hansen, and Kozbur (2016), which introduced
#' the specific implementation used in this package. Heuristically, the penalty parameters are set at
#' a level high enough so that the absolute value of the score for each regressor must be statistically
#' large relative to its standard error in order for the regressors to be selected.
#'
#' @param K Maximum number of iterations.
#' @param penalty Only "lasso" is supported at the present stage.
#' @param phipost Logical. If \code{TRUE}, the plugin coefficient-specific penalty weights are iteratively
#' calculated using estimates from a post-penalty regression. Otherwise, these are calculated using
#' estimates from a penalty regression.
#' @inheritParams penhdfeppml_int
#'
#' @return An object of class \code{elnet} with the elements described in \link[glmnet]{glmnet}, as
#'     well as the following:
#'     \itemize{
#'         \item \code{mu}: a 1 x \code{length(y)} matrix with the final values of the conditional mean \eqn{\mu}.
#'         \item \code{deviance}.
#'         \item \code{bic}: Bayesian Information Criterion.
#'         \item \code{phi}: coefficient-specific penalty weights.
#'         \item \code{x_resid}: matrix of demeaned regressors.
#'         \item \code{z_resid}: vector of demeaned (transformed) dependent variable.
#' }
#'
#' @export
#'
#' @examples
#' # To reduce run time, we keep only countries in the Americas:
#' americas <- countries$iso[countries$region == "Americas"]
#' trade <- trade[(trade$imp %in% americas) & (trade$exp %in% americas), ]
#' # Now generate the needed x, y and fes objects:
#' y <- trade$export
#' x <- data.matrix(trade[, -1:-6])
#' fes <- list(exp_time = interaction(trade$exp, trade$time),
#'             imp_time = interaction(trade$imp, trade$time),
#'             pair     = interaction(trade$exp, trade$imp))
#' # Finally, we try penhdfeppml_cluster_int:
#' reg <- penhdfeppml_cluster_int(y = y, x = x, fes = fes, cluster = fes$pair)
#'
#' @inheritSection hdfeppml_int References

penhdfeppml_cluster_int <- function(y, x, fes, cluster, tol = 1e-8, hdfetol = 1e-4, glmnettol = 1e-12,
                                penalty = "lasso", penweights = NULL, saveX = TRUE, mu = NULL,
                                colcheck_x = TRUE, colcheck_x_fes = TRUE, K = 15, init_z = NULL, post = FALSE,
                                verbose = FALSE, lambda = NULL, phipost=TRUE, gamma_val=NULL) {

  xnames <- colnames(x)
  n <- length(y)
  k <- ncol(x) # BUG? should be defined after colcheck
  nclusters <- nlevels(droplevels(cluster, exclude = if(anyNA(levels(cluster))) NULL else NA))
  x <- data.matrix(x)

  if(is.null(lambda)){
    c <- 1.1
    if(is.null(gamma_val)){gamma_val <- 0.1/log(n)}
    gamma <- gamma_val
    lambda <- c * sqrt(n) * stats::qnorm(1 - gamma / (2 * k))
  }

  if (verbose == TRUE) {
    print("verbose")
  }

  b <- matrix(NA, nrow = ncol(x), ncol = 1)  # fix this later.
  rownames(b) <- colnames(x)
  include_x <- 1:ncol(x)

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

  # number of obs (needed for deviance)
  n <- length(y)

  # estimation algorithm
  crit <- 1
  old_deviance <- 0
  iter <- 0
  while (crit > tol) {
    iter <- iter + 1

    if(iter > 200){message("Plugin Lasso exceeded 200 iterations. Break loop and return last model."); break}

    if (iter == 1) {

      # initialize "mu"
      if (is.null(mu)){
        only_fes <- hdfeppml_int(y, fes=fes, tol = 1e-8, hdfetol = 1e-4, colcheck_x = FALSE, colcheck_x_fes = FALSE, mu = NULL, saveX = TRUE,
                                 init_z = NULL, verbose = FALSE, maxiter = 1000, cluster = NULL, vcv = TRUE)
        mu <- only_fes$mu
      }
      z   <- (y - mu) / mu + log(mu)
      eta <- log(mu)
      last_z <- z
      if (is.null(init_z)) {
        reg_z  <- matrix(z)
      } else{
        reg_z <- init_z
      }
      reg_x  <- x

    } else {
      last_z <- z
      z <- (y - mu) / mu + log(mu)
      reg_z  <- matrix(z - last_z + z_resid)
      reg_x  <- x_resid
      ## colnames(reg_x)   <- colnames(x)
    }

    if(!missing(fes)){
      if(is.null(fes)){
        z_resid <- collapse::fwithin(x=reg_z, g=factor(rep(1,length(reg_z))), w = mu)
        if(!missing(x)){
          x_resid <- collapse::fwithin(x=reg_x,g=factor(rep(1,length(reg_z))), w = mu)
        }
      }else{
        z_resid <-  collapse::fhdwithin(reg_z, fes, w = mu)
        if(!missing(x)){
          x_resid <- collapse::fhdwithin(reg_x, fes, w = mu)
        }
      }
    } else {
      z_resid <- reg_z
      if(!missing(x)){
        x_resid <- reg_x
      }
    }
    # the "cluster_matrix" command computes the variance of the score based on the assumed clustering
    if (iter == 1) {
      e <- mu * z_resid
      phi <- sqrt(diag(cluster_matrix(mu * z_resid, cluster, x_resid)) / n)
    }
    else if (iter<=K) {
      if (phipost == TRUE) {

        x_select <- x_resid[,as.numeric(penreg$beta)!=0]

        if(length(x_select)!=0){
          ppml_temp <- hdfeppml_int(y = y, x = x_select, fes = fes, tol = tol, hdfetol = hdfetol,
                       mu = penreg$mu, colcheck_x = colcheck_x_post, colcheck_x_fes = colcheck_x_fes_post, cluster = cluster)

          mu_post      <- ppml_temp$mu
          x_resid_post <- collapse::fhdwithin(reg_x, fes, w = mu_post)
          residuals_post <- y - mu_post
        }
        else{ # no covariates selected
          residuals_post <-  mu*residuals
          x_resid_post <- x_resid
          mu_post <- mu
        }

        phi <- sqrt(diag(cluster_matrix(residuals_post, cluster, x_resid_post))/n)
      }
      else {
        phi <- sqrt(diag(cluster_matrix(mu*residuals,cluster, x_resid))/n)
      }
    }
    lambda_glmnet <- lambda / sum(y) * sum(phi) / k

    if (verbose == TRUE) {
      print(phi)
    }

# The SCAD option is DISABLED. Since lasso becomes the only option, the if else structure is not needed:
#    if (penalty == "SCAD") {  ## !! currently *very* slow... could be sped up using warm starts.
#      wz_resid <- sqrt(mu) * z_resid
#      wx_resid <- sqrt(mu) * x_resid
#      penreg <- ncvreg::ncvreg(wx_resid, wz_resid, penalty="SCAD", lambda = lambda, penalty.loadings = phi)
#    } else {
      #wx_resid <- cbind(sqrt(mu), wx_resid)
      #penreg <- cdreg(x = x_resid, y = z_resid, weights = mu / sum(mu),lambda = lambda, thresh = 1e-20) ## SLOW!
      if (penalty != "lasso" & iter == 1) {
        warning(penalty, " penalty is not supported. Lasso is used by default.")
      }
      penreg <- glmnet::glmnet(x = x_resid, y = z_resid, weights = mu / sum(mu), lambda = lambda_glmnet,
                               thresh = glmnettol, penalty.factor = phi, standardize = FALSE, family=gaussian(link = "identity"), warm.g=NULL)
#    }

    b[include_x] <- penreg$beta #[-1,]   #does using [,include_x] make a difference?

    if (verbose == TRUE){
      print(penreg$beta)
    }

    residuals <- z_resid - x_resid %*% b[include_x]  #technically this reflects (y-mu)/mu

    mu <- as.numeric(exp(z - residuals))
    mu[which(mu < 1e-190)] <- 1e-190
    mu[mu > 1e190] <- 1e190

    # print("mu")
    # print(length(mu[which(mu == 1e-190)]))

    # calculate deviance
    temp <-  -(y * log(y / mu) - (y - mu))
    temp[which(y == 0)] <- -mu[which(y == 0)]
    #temp[which(mu==Inf)] <- 0

    # print("clust")
    # print(temp[which(is.na(temp))])
    # print(mu[which(is.na(temp))])
    # print(y[which(is.na(temp))])

    deviance <- -2 * sum(temp) / n

    if (deviance < 0) deviance = 0

    delta_deviance <- old_deviance - deviance

    if (!is.na(delta_deviance) & (deviance < 0.1 * delta_deviance)) {
      delta_deviance = deviance
    }

    denom_crit <-  max(c(min(c(deviance, old_deviance)), 0.1))
    crit <-  abs(delta_deviance) / denom_crit

    #print(deviance)

    old_deviance <- deviance
  }

  ## elements to return
  k   <- ncol(matrix(x))
  n   <- length(y)
  select_x <- which(b!=0)

  k   <- length(select_x)
  message(paste("No of variables:", k))
  bic <- deviance + k * log(n) / n
  message(paste("BIC:", bic))
  # k = number of elements in x here
  # BIC would be BIC = deviance + k * ln(n)

  # use ppml estimates instead of lasso estimates
  if (post) {
    x_select <- x_resid[, as.numeric(penreg$beta) != 0]
    if (length(x_select) != 0) {
      ppml_temp <- hdfeppml_int(y = y, x = x_select, fes = fes, tol = tol, hdfetol = hdfetol,
                            mu = penreg$mu, colcheck_x = colcheck_x_post, colcheck_x_fes = colcheck_x_fes_post, cluster = cluster)

      penreg$beta[which(penreg$beta != 0), 1]  <- ppml_temp$coefficients
      b[include_x] <- penreg$beta
      mu    <- ppml_temp$mu
      bic   <- ppml_temp$bic
      deviance <- ppml_temp$deviance
    }
    else{
      print("no covariates selected!")
    }
  }

  penreg[["beta"]] <- b
  penreg[["mu"]]  <-  mu
  penreg[["bic"]] <-  bic
  penreg[["deviance"]] <- deviance
  penreg[["phi"]] <- phi
  penreg[["lambda"]] <- lambda / sum(y)

  if (saveX == TRUE) {
    penreg[["x_resid"]] <- x_resid
    penreg[["z_resid"]] <- z_resid
  }
  return(penreg)
}
