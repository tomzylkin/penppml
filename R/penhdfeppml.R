#' One-Shot Penalized PPML Estimation
#'
#' \code{penhdfeppml} computes a penalized PPML model for a given type of penalty and a given
#' value of the penalty parameter.
#'
#' @param lambda Penalty parameter (a number).
#' @param glmnettol Tolerance parameter to be passed on to glmnet.
#' @param penalty A string. Currently supported: "lasso", "ridge", "SCAD".
#' @param penweights TODO: check what this does (parameter-specific penalties in plugin lasso?).
#' @param post Logical. If \code{TRUE}, estimates a post-penalty model with the selected variables.
#' @param standardize Logical. If \code{TRUE}, x variables are standardized (TODO: check this).
#' @param method TODO: check what this does.
#' @param debug TODO: check what this does.
#' @inheritParams hdfeppml
#'
#' @return A list (TODO: complete this).
#' @export
#'
#' @examples # TODO: add examples here.

penhdfeppml <- function(y, x, fes, lambda, tol = 1e-8, hdfetol = 1e-4, glmnettol = 1e-12,
                        penalty = "lasso", penweights = NULL,
                        selectobs = rep(TRUE, length(y)), saveX = TRUE, mu = NULL, colcheck = TRUE,
                        init_z = NULL, post = FALSE, verbose = FALSE, standardize = TRUE,
                        method = "placeholder", cluster = NULL, debug = FALSE) {

  # We'll subset using selectobs up front (no need to keep calling selectobs later on)
  if (!is.null(selectobs)) {
    y <- y[selectobs]
    x <- x[selectobs, ] # Subsetting x. This works even if x is a vector: we coerced it to matrix in l.56.
    # Subsetting fes (we're using a for loop because we're assuming they're in list form):
    for (i in seq_along(fes)) {
      fes[[i]] <- fes[[i]][selectobs]
    }
    # Important: we need to subset clusters too (if used):
    if (!is.null(cluster)) cluster <- cluster[selectobs]
  }

  # implements plugin method; calls penhdfeppml_cluster subcommand (careful: selectobs is NULL, because
  # we already filtered y, x, fes and cluster):
  if (method == "iterative") {
    penreg <- penhdfeppml_cluster(y = y, x = x, fes = fes, cluster = cluster, tol = tol,
                                  hdfetol = hdfetol, glmnettol = glmnettol, penalty = penalty,
                                  penweights = penweights, selectobs = NULL, saveX = saveX,
                                  mu = mu, colcheck = colcheck, K = 15, init_z = init_z, post = FALSE,
                                  verbose = verbose, lambda = NULL)
  }
  else {

    # if "iterative" option not enabled, do the following
    b <- matrix(NA, nrow = ncol(x), ncol = 1)  # fix this later.
    rownames(b) <- colnames(x)
    include_x <- 1:ncol(x)

    if (is.null(penweights)) {
      penweights <- rep(1, length(include_x))   #note this is the default used by glmnet. Using "NULL" actually gives you incorrect weights.
    }

    if (length(penweights) > length(include_x)){
      print("penweights needs to be same length as number of x variables")
      stop()
    }

    # collinearity check
    if (colcheck == TRUE) {
      include_x <- collinearity_check(y, x, fes, 1e-6)
      x <- x[, include_x]
      if (!is.null(penweights)) {
        penweights = penweights[include_x]
      }
    }

    # number of obs (needed for deviance)
    n <- length(y)

    # estimation algorithm (this implements a version of the algorithm from p. 110 of CGZ SJ 2020 but with penalized WLS in place of the WLS step)
    crit <- 1
    old_deviance <-0
    iter <-0
    while (crit > tol) {
      iter <- iter + 1

      if (iter == 1) {

        # initilize "mu"
        if (is.null(mu)) mu  <- (y + mean(y))/2
        z   <- (y - mu)/mu + log(mu)
        eta <- log(mu)
        last_z <- z
        if (is.null(init_z)) {
          reg_z  <- matrix(z)
        } else {
          reg_z <- init_z
        }
        reg_x  <- x

      } else {
        last_z <- z
        z <- (y - mu)/mu + log(mu)
        reg_z  <- matrix(z - last_z + z_resid)
        reg_x  <- x_resid
        ## colnames(reg_x)   <- colnames(x)
      }

      # HDFE estimation works with the residuals of z and x after purging them of the FEs (see CGZ Stata Journal 2020)
      z_resid <- lfe::demeanlist(reg_z, fes, weights = sqrt(mu), eps = hdfetol)
      x_resid <- lfe::demeanlist(reg_x, fes, weights = sqrt(mu), eps = hdfetol)

      if (is.null(penweights)) {
        #penreg <- glmnet::glmnet(x = x_resid, y = z_resid, weights = mu/sum(mu), lambda = lambda, thresh = glmnettol, standardize = standardize)
      }
      else{
        if(debug) {
          print(sum(penweights))
        }
        if(debug) {
          print(penweights)
        }
        penweights = penweights * length(include_x) / sum(penweights)
        if(debug) {
          print(sum(penweights))
        }
      }

      if (penalty == "SCAD") {  ## !! currently *very* slow... can be sped up using warm starts??
        wz_resid <- sqrt(mu) * z_resid
        wx_resid <- sqrt(mu) * x_resid
        penreg <- ncvreg::ncvreg(wx_resid, wz_resid, penalty = "SCAD", lambda = lambda)  # add penalty weights

      }  else if (penalty == "ridge") {
        penreg <- fastridge(x = x_resid, y = z_resid, weights = mu/sum(mu), lambda = n * lambda,
                            standardize = standardize) # ,penalty.factor=penweights
      } else {

        # Lasso is the default
        if (debug) {
          penreg <- glmnet::glmnet(x = x_resid, y = z_resid, weights = mu/sum(mu), lambda = lambda,
                                   thresh = glmnettol, standardize = standardize)
          print((penreg$beta))
          print(penweights)
        }
        if (debug) {
          penreg <- glmnet::glmnet(x = x_resid, y = z_resid, weights = mu/sum(mu), lambda = lambda,
                                   thresh = glmnettol, standardize = standardize)
          print((penreg$beta))
          print(penweights)
        }
        penreg <- glmnet::glmnet(x = x_resid, y = z_resid, weights = mu/sum(mu), lambda = lambda,
                                 thresh = glmnettol, penalty.factor = penweights, standardize = standardize)
        if (debug) {
          print((penreg$beta))
          stop()
        }
      }

      b[include_x] <- penreg$beta #[-1,]   #does using [,include_x] make a difference?

      residuals <- z_resid - x_resid %*% b[include_x]

      mu <- as.numeric(exp(z - residuals))

      # calculate deviance
      temp <-  -(y * log(y/mu) - (y-mu))
      temp[which(y == 0)] <- -mu[which(y == 0)]
      if (!missing(selectobs)){
        temp[which(!selectobs)] <-0
      }

      deviance <- -2 * sum(temp)/n

      if(deviance<0) deviance = 0

      delta_deviance <- old_deviance - deviance

      if (!is.na(delta_deviance) & (deviance < 0.1 * delta_deviance)) {
        delta_deviance = deviance
      }

      denom_crit = max(c( min(c(deviance, old_deviance)) , 0.1 ))
      crit = abs(delta_deviance) / denom_crit

      #print(deviance)
      if (verbose == TRUE) {
        #print(deviance)
        print(crit)
      }

      old_deviance <- deviance
    }

    ## elements to return
    k   <- ncol(matrix(x))
    n   <- length(y)
    select_x <- which(b != 0)

    k   <- length(select_x) #ncol(matrix(x[,select_x]))
    bic <- deviance + k * log(n)/n
    # k = number of elements in x here
    # BIC would be BIC = deviance + k * ln(n)


    # if "post" enabled use post-lasso ppml estimates instead of lasso estimates
    if (post) {
      x_select <- x_resid[, as.numeric(penreg$beta) != 0]
      if (length(x_select) != 0){
        ppml_temp <- hdfeppml(y = y, x = x_select, fes = fes, tol = tol, hdfetol = hdfetol,
                              mu = penreg$mu, colcheck = FALSE)

        penreg$pencoefs <- penreg$beta
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

    # return results
    penreg[["beta"]] <- b
    penreg[["mu"]]  <-  mu
    penreg[["bic"]] <-  bic
    penreg[["deviance"]] <- deviance
    if (saveX == TRUE) {
      penreg[["x_resid"]] <- x_resid
      penreg[["z_resid"]] <- z_resid
    }
  }
  return(penreg)
}
