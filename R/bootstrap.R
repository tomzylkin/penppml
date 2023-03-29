#' Bootstrap Lasso Implementation (in development)
#'
#' This function performs standard plugin lasso PPML estimation for \code{bootreps} samples drawn again with
#' replacement and reports
#' those regressors selected in at least a certain fraction of the bootstrap repetitions.
#'
#' This function enables users to implement the "bootstrap" step in the procedure described in
#' Breinlich, Corradi, Rocha, Ruta, Santos Silva and Zylkin (2020). To do this, Plugin Lasso is run B times.
#' The function can also perform a post-selection estimation.
#'
#' @param data A data frame containing all relevant variables.
#' @param dep A string with the names of the independent variables or their column numbers.
#' @param indep A vector with the names or column numbers of the regressors. If left unspecified,
#'              all remaining variables (excluding fixed effects) are included in the regressor matrix.
#' @param cluster_id A string denoting the cluster-id with which to perform
#'  cluster bootstrap.
#' @param selectobs Optional. A vector indicating which observations to use (either a logical vector
#'     or a numeric vector with row numbers, as usual when subsetting in R).
#' @param bootreps Number of bootstrap repetitions.
#' @param boot_threshold Minimal threshold. If a variable is selected in at least this
#' fraction of times, it is reported at the end of the iterations.
#' @param colcheck_x Logical. If \code{TRUE}, this checks collinearity between the independent variables and drops the
#' collinear variables.
#' @param colcheck_x_fes Logical. If \code{TRUE}, this checks whether the independent variables are perfectly explained
#' by the fixed effects drops those that are perfectly explained.
#' @param phipost Logical. If \code{TRUE}, the plugin coefficient-specific penalty weights are iteratively
#' calculated using estimates from a post-penalty regression. Otherwise,
#' these are calculated using estimates from a penalty regression.
#'@inheritParams mlfitppml
#'@inheritParams penhdfeppml_int
#'@inheritParams hdfeppml_int
#'
#' @return A matrix with coefficient estimates for all dependent variables.
#' @export
#'
#' @examples
#' \dontrun{bs1 <- bootstrap(data=trade3, dep="export",
#'                  cluster_id="clus",
#'                  fixed=list(c("exp", "time"),
#'                  c("imp", "time"), c("exp", "imp")),
#'                  indep=7:22, bootreps=10, colcheck_x = TRUE,
#'                  colcheck_x_fes = TRUE,
#'                  boot_threshold = 0.01,
#'                  post=TRUE, gamma_val=0.01, verbose=FALSE)}
#'
#' @inheritSection hdfeppml_int References

bootstrap <- function(data, dep, indep = NULL, cluster_id=NULL, fixed=NULL, selectobs = NULL, bootreps=250, boot_threshold = 0.01, colcheck_x=FALSE,
                      colcheck_x_fes=FALSE, post=FALSE, gamma_val = NULL,
                      verbose=FALSE, tol=1e-6, hdfetol=1e-2, penweights=NULL, maxiter=1000, phipost=TRUE){
  cluster <- one <- NULL
  # First we do the data handling with genmodel:
  model <- genmodel(data = data, dep = dep, indep = indep, selectobs = selectobs)

  indep_names <- colnames(data[,indep])

  # Now we create the result matrices:
  bootstrap_results <- list()
  uniq_clus <- levels(factor(data[,cluster_id])) # Unique clusters, as characters
  save_betas     <- matrix(nrow = length(indep), ncol = bootreps, dimnames = list(colnames(data[indep]),1:bootreps))
  save_betas_pre <- matrix(nrow = length(indep), ncol = bootreps, dimnames = list(colnames(data[indep]),1:bootreps))
  save_phis      <- matrix(nrow = length(indep), ncol = bootreps, dimnames = list(colnames(data[indep]),1:bootreps))
  is_included    <- matrix(nrow = length(indep), ncol = bootreps, dimnames = list(colnames(data[indep]),1:bootreps))
  draws     <- matrix(nrow = length(uniq_clus), ncol = bootreps)

  for (b in 1:bootreps) { # draw matrix has resampled cluster IDs in it
    draws[,b] <- sort(sample(uniq_clus,replace=TRUE))
  }

  # Run bootstrap repetitions
  for (b in 1:bootreps) {

    # Close sink file if it wasn't closed in the last bootstrap repetition
    if(b>1){
      if(all(is.na(is_included[,b-1]))){
        sink()
      }
    }

    tryCatch({
    draw  <- draws[,b] # Take one draw
    draw <- data.frame(cbind(draw,1))
    colnames(draw) <- c("cluster","one")
    draw <- data.frame(draw) %>% dplyr::group_by(cluster) %>% dplyr::mutate(N=cumsum(one))
    draw <- draw[,-2]
    colnames(draw)[1] <- cluster_id # First column contains cluster IDs
    boot_data  <- merge(data, draw, cluster_id, all.x=FALSE, all.y=TRUE) # Merge (cluster, data) and draw, keep only those drawn, i.e. those that appear in draw
    boot_index <-  rep(row.names(boot_data), boot_data$N) # Repeat row name the times indicated in draw
    boot_data  <- boot_data[boot_index, ] # Get the row names indicated like this
    boot_data$rep <- (as.numeric(rownames(boot_data)) %% 1)*10+1 # This creates an ID for
    boot_rep  <- factor(boot_data$rep)

    #boot_ID = boot_data$id
    #boot_ID[is.na(WB_TRADE_DATA$alt_id)] <- 0
    #boot_data$pair <- interaction(boot_exp, boot_imp)
    #boot_data <- within(boot_data, {alt_id2 = ave(alt_id,pair,FUN=max)} )
    #boot_data$alt_id2[boot_data$alt_id2==0] <- boot_data$pair[boot_data$alt_id2==0]

    boot_clus    <- boot_data[,cluster_id]
    boot_clus    <- interaction(boot_clus, boot_rep)
    boot_data$clus2 <- boot_clus

    boot_x        <- boot_data[,indep_names]

    #boot_mu       <- boot_data$mu
    boot_rta      <- ((rowSums(boot_x ))>0)*1
    boot_data$rta <- boot_rta
    boot_rta_only <- hdfeppml(dep=dep, indep="rta", fixed = fixed,
                              tol=tol,hdfetol=hdfetol,cluster="clus2", colcheck_x=colcheck_x,
                              colcheck_x_fes=colcheck_x_fes, data=boot_data, verbose=verbose, maxiter=maxiter)
    boot_mu <- boot_rta_only$mu

    # This runs bootstrap repetition of plugin
    sink(file="NULL")
    plugin_boot <- mlfitppml(dep=dep,indep=indep_names,fixed = fixed, tol=tol, hdfetol=hdfetol,
                                            cluster="clus2",method="plugin", colcheck_x=colcheck_x,
                                            colcheck_x_fes=colcheck_x_fes, mu=boot_mu, data=boot_data,
                                            gamma_val = gamma_val, verbose=verbose, phipost=phipost, penweights=penweights)
    sink()
    # add option for orig.
    save_betas[rownames(plugin_boot$beta),b]     <- plugin_boot$beta
    save_betas_pre[rownames(plugin_boot$beta),b] <- plugin_boot$beta_pre
    save_phis[rownames(plugin_boot$beta),b]      <- plugin_boot$phi
    is_included[rownames(plugin_boot$beta),b]    <- t(colSums(boot_data[,rownames(plugin_boot$beta)]) !=0)
    is_included[is.na(is_included[,b]),b] <- FALSE
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }

  selected_vars <- names(which(rowSums(save_betas!=0, na.rm=T)/bootreps>boot_threshold))

  if(post==TRUE){
  if(length(selected_vars)>0){
    boot_post <- hdfeppml(dep=dep, indep=selected_vars, fixed = fixed,
                       tol=1e-6,hdfetol=1e-2,cluster="clus2", colcheck_x=colcheck_x, colcheck_x_fes=colcheck_x_fes, data=boot_data)
  } else {
    message("No variable selected by Bootstrap Lasso and therefore cannot compute Post-Bootstrap Lasso.")
  }
  }

  bootstrap_results[["betas"]] <- save_betas
  bootstrap_results[["betas_pre"]] <- save_betas_pre
  bootstrap_results[["phis"]] <- save_phis
  bootstrap_results[["is_included"]] <- is_included
  bootstrap_results[["selected"]] <- selected_vars
 if(post==TRUE & length(selected_vars)>0){
    bootstrap_results[["post"]] <- boot_post$coefficients
   }

  return(bootstrap_results)
}
