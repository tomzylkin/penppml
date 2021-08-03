#' Plugin Lasso Estimation
#'
#' Estimates coefficient-specific penalty weights that account for heteroskedasticity.
#' Called by \code{mlfitppml} and \code{penhdfeppml}.
#'
#' @param K TODO: check what this does.
#' @inheritParams penhdfeppml
#'
#' @return TODO: check this.
#' @export
#'
#' @examples # TODO: add examples here.

penhdfeppml_cluster <- function(y, x, fes, cluster, tol = 1e-8, hdfetol = 1e-4, glmnettol = 1e-12,
                                penalty = "lasso", penweights = NULL, selectobs = NULL, saveX = TRUE,
                                mu = NULL, colcheck = TRUE, K = 15, init_z = NULL, post = FALSE,
                                verbose = FALSE, lambda = NULL) {

  # allow for more penalty functions

  n <- length(y)
  k <- ncol(x) # BUG? should be defined after colcheck
  nclusters <- nlevels(droplevels(cluster, exclude = if(anyNA(levels(cluster))) NULL else NA))
  x <- data.matrix(x)

  if(is.null(lambda)){
    c <- 1.1
    gamma <- .1 / log(n)
    lambda <- c * sqrt(n) * stats::qnorm(1 - gamma / (2 * k))
  }

  # We'll subset using selectobs up front (no need to keep calling selectobs later on)
  if (!is.null(selectobs)) {
    y <- y[selectobs]
    x <- x[selectobs, ] # Subsetting x. This works even if x is a vector: we coerced it to matrix in l.24.
    # Subsetting fes (we're using a for loop because we're assuming they're in list form):
    for (i in seq_along(fes)) {
      fes[[i]] <- fes[[i]][selectobs]
    }
    # Important: we need to subset clusters too:
    cluster <- cluster[selectobs]
  }

  if(verbose == TRUE){
    print("verbose")
  }

  b <- matrix(NA, nrow = ncol(x), ncol = 1)  # fix this later.
  rownames(b) <- colnames(x)
  include_x <- 1:ncol(x)

  if (colcheck == TRUE){
    include_x <- collinearity_check(y, x, fes, 1e-6)
    x <- x[, include_x]
  }

  # number of obs (needed for deviance)
  n   <- length(y)

  # estimation algorithm
  crit <- 1
  old_deviance <- 0
  iter <- 0
  while (crit > tol) {
    iter <- iter + 1

    if (iter == 1) {

      # initilize "mu"
      if (is.null(mu)) mu  <- (y + mean(y)) / 2
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
    print(iter)

    z_resid <- lfe::demeanlist(reg_z, fes, weights = sqrt(mu), eps = hdfetol)
    x_resid <- lfe::demeanlist(reg_x, fes, weights = sqrt(mu), eps = hdfetol)

    # the "cluster_matrix" command computes the variance of the score based on the assumed clustering
    if (iter == 1) {
      e <- mu * z_resid
      phi <- sqrt(diag(cluster_matrix(mu * z_resid, cluster, x_resid)) / n)
    }
    else if (iter<=K) {
      phi <- sqrt(diag(cluster_matrix(mu * residuals, cluster, x_resid)) / n) # alternatively can use y-mu
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
                               thresh = glmnettol, penalty.factor = phi, standardize = FALSE)
#    }

    b[include_x] <- penreg$beta #[-1,]   #does using [,include_x] make a difference?

    if (verbose == TRUE){
      print(penreg$beta)
    }

    residuals <- z_resid - x_resid %*% b[include_x]  #technically this reflects (y-mu)/mu

    mu <- as.numeric(exp(z - residuals))

    # calculate deviance
    temp <-  -(y * log(y / mu) - (y - mu))
    temp[which(y == 0)] <- -mu[which(y == 0)]
    if (!missing(selectobs)){
      temp[which(!selectobs)] <- 0
    }

    deviance <- -2 * sum(temp) / n

    if(deviance < 0) deviance = 0

    delta_deviance <- old_deviance - deviance

    if (!is.na(delta_deviance) & (deviance < 0.1 * delta_deviance)) {
      delta_deviance = deviance
    }

    denom_crit = max(c(min(c(deviance, old_deviance)), 0.1 ))
    crit = abs(delta_deviance) / denom_crit

    #print(deviance)

    old_deviance <- deviance
  }

  ## elements to return
  k   <- ncol(matrix(x))
  n   <- length(y)
  select_x <- which(b!=0)

  k   <- length(select_x)
  print(k)
  bic <- deviance + k * log(n) / n
  # k = number of elements in x here
  # BIC would be BIC = deviance + k * ln(n)

  # use ppml estimates instead of lasso estimates
  if (post) {
    x_select <- x_resid[, as.numeric(penreg$beta) != 0]
    if(length(x_select) != 0){
      ppml_temp <- hdfeppml(y = y, x = x_select, fes = fes, tol = tol, hdfetol = hdfetol,
                            mu = penreg$mu, colcheck = FALSE, cluster = cluster)

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
