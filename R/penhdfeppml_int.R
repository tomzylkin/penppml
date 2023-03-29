#' One-Shot Penalized PPML Estimation with HDFE
#'
#' \code{penhdfeppml_int} is the internal algorithm called by \code{penhdfeppml} to fit a penalized PPML
#' regression for a given type of penalty and a given value of the penalty parameter. It takes a vector
#' with the dependent variable, a regressor matrix and a set of fixed effects (in list form: each element
#' in the list should be a separate HDFE). The penalty can be either lasso or ridge, and the plugin
#' method can be enabled via the \code{method} argument.
#'
#' More formally, \code{penhdfeppml_int} performs iteratively re-weighted least squares (IRLS) on a
#' transformed model, as described in Breinlich, Corradi, Rocha, Ruta, Santos Silva and Zylkin (2020).
#' In each iteration, the function calculates the transformed dependent variable, partials out the fixed
#' effects (calling \code{collapse::fhdwithin}) and then and then calls \code{glmnet} if the selected
#' penalty is lasso (the default). If the user selects ridge, the analytical solution is instead
#' computed directly using fast C++ implementation.
#'
#' For information on the plugin lasso method, see \link{penhdfeppml_cluster_int}.
#'
#' @param lambda Penalty parameter (a number).
#' @param glmnettol Tolerance parameter to be passed on to \code{glmnet}.
#' @param penalty A string indicating the penalty type. Currently supported: "lasso" and "ridge".
#' @param penweights Optional: a vector of coefficient-specific penalties to use in plugin lasso when
#'     \code{method == "plugin"}.
#' @param post Logical. If \code{TRUE}, estimates a post-penalty regression with the selected variables.
#' @param standardize Logical. If \code{TRUE}, x variables are standardized before estimation.
#' @param method The user can set this equal to "plugin" to perform the plugin algorithm with
#'     coefficient-specific penalty weights (see details). Otherwise, a single global penalty is used.
#' @param debug Logical. If \code{TRUE}, this helps with debugging penalty weights by printing output
#'    of the first iteration to the console and stopping the estimation algorithm.
#' @param colcheck_x Logical. If \code{TRUE}, this checks collinearity between the independent variables and drops the
#' collinear variables.
#' @param colcheck_x_fes Logical. If \code{TRUE}, this checks whether the independent variables are perfectly explained
#' by the fixed effects drops those that are perfectly explained.
#' @param gamma_val Numerical value that determines the regularization threshold as defined in Belloni, Chernozhukov, Hansen, and Kozbur (2016). NULL default sets parameter to 0.1/log(n).
#' @param phipost Logical. If \code{TRUE}, the plugin coefficient-specific penalty weights are iteratively
#' calculated using estimates from a post-penalty regression when \code{method == "plugin"}. Otherwise,
#' these are calculated using estimates from a penalty regression.
#' @inheritParams hdfeppml_int
#'
#' @return If \code{method == "lasso"} (the default), an object of class \code{elnet} with the elements
#'     described in \link[glmnet]{glmnet}, as well as:
#'     \itemize{
#'         \item \code{mu}: a 1 x \code{length(y)} matrix with the final values of the conditional mean \eqn{\mu}.
#'         \item \code{deviance}.
#'         \item \code{bic}: Bayesian Information Criterion.
#'         \item \code{phi}: coefficient-specific penalty weights (only if \code{method == "plugin"}.
#'         \item \code{x_resid}: matrix of demeaned regressors.
#'         \item \code{z_resid}: vector of demeaned (transformed) dependent variable.
#' }
#'     If \code{method == "ridge"}, a list with the following elements:
#' \itemize{
#'   \item \code{beta}: a 1 x \code{ncol(x)} matrix with coefficient (beta) estimates.
#'   \item \code{mu}: a 1 x \code{length(y)} matrix with the final values of the conditional mean \eqn{\mu}.
#'   \item \code{deviance}.
#'   \item \code{bic}: Bayesian Information Criterion.
#'   \item \code{x_resid}: matrix of demeaned regressors.
#'   \item \code{z_resid}: vector of demeaned (transformed) dependent variable.
#' }
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
#' # Finally, we try penhdfeppml_int with a lasso penalty (the default):
#' reg <- penhdfeppml_int(y = y, x = x, fes = fes, lambda = 0.1)
#'
#' # We can also try ridge:
#' \donttest{reg <- penhdfeppml_int(y = y, x = x, fes = fes, lambda = 0.1, penalty = "ridge")}
#'
#' @inheritSection hdfeppml_int References
#' @importFrom stats gaussian
#' @importFrom stats var
#' @importFrom utils head
#' @importFrom devtools load_all

penhdfeppml_int <- function(y, x, fes, lambda, tol = 1e-8, hdfetol = 1e-4, glmnettol = 1e-12,
                            penalty = "lasso", penweights = NULL, saveX = TRUE, mu = NULL, colcheck_x = FALSE, colcheck_x_fes = TRUE,
                            init_z = NULL, post = FALSE, verbose = FALSE, phipost = TRUE, standardize = TRUE,
                            method = "placeholder", cluster = NULL, debug = FALSE, gamma_val=NULL) {
  xnames <- colnames(x)
  old_x <- x
  old_y <- y
  old_fes <- fes
  # implements plugin method; calls penhdfeppml_cluster_int subcommand
  if (method == "plugin") {
    penreg <- penhdfeppml_cluster_int(y = y, x = x, fes = fes, cluster = cluster, tol = tol,
                                      hdfetol = hdfetol, glmnettol = glmnettol, penalty = penalty,
                                      penweights = penweights, saveX = saveX,
                                      mu = mu, colcheck_x = colcheck_x, colcheck_x_fes = colcheck_x_fes, K = 15, init_z = init_z, post = FALSE,
                                      verbose = verbose, phipost=phipost, lambda = NULL, gamma_val=gamma_val)
  }
  else {

    # if "plugin" option not enabled, do the following
    b <- matrix(NA, nrow = ncol(x), ncol = 1)  # fix this later.
    rownames(b) <- colnames(x)
    include_x <- 1:ncol(x)

   # if (is.null(penweights)) {
   #    penweights <- rep(1, length(include_x))   #note this is the default used by glmnet. Using "NULL" actually gives you incorrect weights.
   #  } #N: moved this below colchecks, because include_x gets updated

    if (length(penweights) > length(include_x)){
      print("penweights needs to be same length as number of x variables")
      stop()
    }

    # collinearity check
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
    if (colcheck_x==FALSE & colcheck_x_fes == TRUE){
      include_x <- collinearity_check(y,x,fes,1e-6, colcheck_x=colcheck_x, colcheck_x_fes=colcheck_x_fes)
      x <- x[,include_x]
      colnames(x) <- xnames[include_x]
      xnames <- xnames[include_x]
      colcheck_x_post = FALSE
      colcheck_x_fes_post = TRUE
    }

    if (is.null(penweights)) {
      penweights <- rep(1, length(include_x))   #note this is the default used by glmnet. Using "NULL" actually gives you incorrect weights.
    }

    # number of obs (needed for deviance)
    n <- length(y)
    # estimation algorithm (this implements a version of the algorithm from p. 110 of CGZ SJ 2020 but with penalized WLS in place of the WLS step)
    crit <- 1
    old_deviance <-0
    iter <-0
    while (crit > tol) {
      iter <- iter + 1

      if(iter > 50){message("Lasso exceeded 50 iterations. Break loop and return last model."); break}

      if (iter == 1) {
        # initilize "mu"
        if (is.null(mu)){
          only_fes <- hdfeppml_int(y, fes=fes, tol = 1e-8, hdfetol = 1e-4, colcheck_x = FALSE, colcheck_x_fes = FALSE, mu = NULL, saveX = TRUE,
                                   init_z = NULL, verbose = FALSE, maxiter = 1000, cluster = NULL, vcv = TRUE)
          mu <- only_fes$mu
          #      mu <- mu[mu>0]
          # print(length(mu))
        }
        z   <- (y - mu)/mu + log(mu)
        #     z <- z[mu>0]
        eta <- log(mu)
        last_z <- z
        if (is.null(init_z)) {
          reg_z  <- matrix(z)
        } else {
          reg_z <- init_z
        }
        reg_x  <- x[which(mu>0),]
        n <- length(z)

      } else {
        last_z <- z
        n <- length(mu)
        # print(n)
        #        y <- y[mu>0]
        z <- (y - mu)/mu + log(mu)
        #z[which(mu==1e-190)] <- 1e-190
        #z[which(z==Inf)] <- 1
        reg_z  <- matrix(z - last_z[which(mu>0)] + z_resid[which(mu>0)])
        #print("reg_z")
        #print(sum(is.na(reg_z)))
        reg_x  <- x_resid[which(mu>0),]
        ## colnames(reg_x)   <- colnames(x)
      }
      #      fes <- lapply(fes, "[", which(mu>0))
      # print(length(reg_z))
      # print(dim(reg_x)[1])
      # print(length(mu))
      # print(length(fes[[1]]))
      # print(length(fes[[2]]))
      # print(length(fes[[3]]))
      # HDFE estimation works with the residuals of z and x after purging them of the FEs (see CGZ Stata Journal 2020)
      if(!missing(fes)){
        #   print("fes not missing")
        if(is.null(fes)){
          z_resid <- collapse::fwithin(x=reg_z, g=factor(rep(1,length(reg_z))), w = mu)
          if(!missing(x)){
            x_resid <- collapse::fwithin(x=reg_x,g=factor(rep(1,length(reg_z))), w = mu)
          }
        }else{
          # print(length(reg_z))
          # print(length(fes[[1]]))
          # print(length(mu))
          z_resid <-  collapse::fhdwithin(reg_z, fes, w = mu)
          if(!missing(x)){
            x_resid <- collapse::fhdwithin(reg_x, fes, w = mu)
          }
          # print("postfhd1")
        }
      } else {
        print("fes missing")
        z_resid <- reg_z
        if(!missing(x)){
          x_resid <- reg_x
        }
      }

      # if(sd(z_resid)==Inf){
      #   print("z_resid")
      #   print(head(z_resid))
      #   print(sum(is.na(z_resid)))
      #   print(mean(z_resid))
      #   print(sd(z_resid))
      # }
      # if(length(reg_z)!=length(z_resid)){
      #   mu_t <- mu
      #   reg_z_t <- reg_z
      #   temp1 <- (collapse::fhdwithin(reg_z, fes, w = mu, eps = hdfetol))
      #   temp2 <- (lfe::demeanlist(reg_z, fes, weights = sqrt(mu), eps = hdfetol))
      #   print("check")
      #   print(all.equal(temp1,temp2))
      #   print(sum(is.na(temp2)))
      # }

      #       print(length(z))
      #       print(length(reg_z))
      #print("mu")
      #print(length(mu))
      #print(reg_z)
      # print(sum(mu==0))
      # print("data")
      #  print(length(z_resid))
      #  print(dim(reg_x))
      # print(length(z_resid))
      # print(dim(x_resid))

      if (is.null(penweights)) {
        print("null penweights")
        #penreg <- glmnet(x = x_resid, y = z_resid, weights = mu/sum(mu), lambda = lambda, thresh = glmnettol, standardize = standardize, family=gaussian(link = "identity"), warm.g=NULL)
      }
      else{
        if(debug) {
          print("sum penweights")
          print(sum(penweights))
        }
        if(debug) {
          print("penweights")
          print(penweights)
        }
        penweights = penweights * length(include_x) / sum(penweights)
        if(debug) {
          print(sum(penweights))
        }
      }

      # The SCAD option is DISABLED:
      #      if (penalty == "SCAD") {  ## !! currently *very* slow... can be sped up using warm starts??
      #        wz_resid <- sqrt(mu) * z_resid
      #        wx_resid <- sqrt(mu) * x_resid
      #        penreg <- ncvreg::ncvreg(wx_resid, wz_resid, penalty = "SCAD", lambda = lambda)  # add penalty weights
      #
      #      }  else if (penalty == "ridge") {
      if (penalty == "ridge") {
        penreg <- fastridge(x = x_resid, y = z_resid, weights = mu/sum(mu), lambda = n * lambda,
                            standardize = standardize) # ,penalty.factor=penweights
      } else {

        # Lasso is the default, so a
        if (penalty != "lasso" & iter == 1) {
          warning(penalty, " penalty is not supported. Lasso is used by default.")
        }
        if (debug) {
          penreg <- glmnet::glmnet(x = x_resid, y = z_resid, weights = mu/sum(mu), lambda = lambda,
                                thresh = glmnettol, standardize = standardize, family=gaussian(link = "identity"))
          print((penreg$beta))
          print(penweights)
        }
        if (debug) {
          penreg <- glmnet::glmnet(x = x_resid, y = z_resid, weights = mu/sum(mu), lambda = lambda,
                                thresh = glmnettol, standardize = standardize, family=gaussian(link = "identity"))
          print((penreg$beta))
          print(penweights)
        }
        penreg <- glmnet::glmnet(x = x_resid, y = z_resid, weights = mu/sum(mu), lambda = lambda,
                              thresh = glmnettol, penalty.factor = penweights, standardize = standardize, family=gaussian(link = "identity"))
        if (debug) {
          print((penreg$beta))
          stop()
        }
      }

      b[include_x] <- penreg$beta #[-1,]   #does using [,include_x] make a difference?

      residuals <- z_resid - x_resid %*% b[include_x]

      mu <- as.numeric(exp(z - residuals))
      #print("smaller zero")
      #     print(length(which(mu <= 0)))
      #      mu <- mu[which(mu > 0)]
      mu[which(mu < 1e-190)] <- 1e-190
      mu[mu > 1e190] <- 1e190
      #y <- y[which(mu > 0)]
      # print("mu")
      # print(length(mu[which(mu == 1e-190)]))
      # print(length(mu[which(mu == 1e16)]))

      # calculate deviance
      temp <-  -(y * log(y/mu) - (y-mu)) # Problem sits here: - Inf + Inf = NaN
      temp[which(y == 0)] <- -mu[which(y == 0)]
      #temp[which(mu==Inf)] <- 0
      # print("temp")
      # print(y[which(is.na(temp))])
      # print(mu[which(is.na(temp))])

      deviance <- -2 * sum(temp)/n
      # print("pen")
      # print(temp[which(is.na(temp))])
      # print(mu[which(is.na(temp))])
      # print(y[which(is.na(temp))])
      if(deviance<0) deviance = 0

      delta_deviance <- old_deviance - deviance

      if (!is.na(delta_deviance) & (deviance < 0.1 * delta_deviance)) {
        delta_deviance = deviance
      }

      denom_crit = max(c( min(c(deviance, old_deviance)) , 0.1 ))
      crit = abs(delta_deviance) / denom_crit

      if (verbose == TRUE) {
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
        ppml_temp <- hdfeppml_int(y = y, x = x_select, fes = fes, tol = tol, hdfetol = hdfetol,
                                  mu = penreg$mu, colcheck_x = colcheck_x_post, colcheck_x_fes = colcheck_x_fes_post)
        penreg$pencoefs <- penreg$beta
        penreg$beta[which(as.logical(penreg$beta != 0)), 1]  <- ppml_temp$coefficients
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
