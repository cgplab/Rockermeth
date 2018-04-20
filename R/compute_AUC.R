#' Compute Area Under Curve for a matrix of samples
#'
#' This function computes the Area Under Curve that is then used to define the
#' segregation between tumor and control samples accordingly to their
#' methylation (beta) values.
#'
#' @param tumor_table A matrix of beta-values (percentage) from tumor samples.
#' @param control_table A matrix of beta-values (percentage) from normal/control samples.
#' @param ncores Number of parallel processes to use for parallel computing
#' @param na_threshold Fraction of NAs (considered independently in tumor and
#' control samples) above which a site will not be selected (default=0)
#' @return A vector of AUC scores
#' @export
compute_AUC <- function(tumor_table, control_table, ncores = 1, na_threshold = 0) {
  ncores <- as.integer(ncores)
  max_cores <- parallel::detectCores()
  if (ncores > max_cores){
    stop(sprintf("Selected %i cores but system has %i cores.",
                 ncores, max_cores))
  }

  beta_table <- as.matrix(cbind(tumor_table, control_table))
  sample_state <- c(rep(T, ncol(tumor_table)), rep(F, ncol(control_table)))
  if (all(beta_table >= 0, na.rm = T) && all(beta_table <= 1, na.rm = T) &&
    is.double(beta_table[1,1])) {
    beta_table <- beta_table*100
    storage.mode(beta_table) <- "integer"
  } else {
    stop("tumor_table and control_table must have fraction values")
  }

  na_threshold <- as.numeric(na_threshold)
  stopifnot(na_threshold >= 0 || na_threshold < 1)

  cl <- parallel::makeCluster(nclust)
  auc <- parallel::parApply(cl, beta_table, 1, single_AUC,
    state = sample_state, na_threshold = na_threshold)
  parallel::stopCluster(cl)
  return(auc)
}

#' Compute AUC a single vector
#'
#' Return NA if NA samples are more than threshold
#'
#' @param x integer vector (range 1-100)
#' @param state logical vector
#' @param na_threshold numeric
#' @keywords internal
single_AUC <- function(x, state, na_threshold) {
  all_NAs <- all(is.na(x))
  t_NA_frac <- sum(is.na(x[state])) / length(x[state]) > na_threshold
  c_NA_frac <- sum(is.na(x[!state])) / length(x[!state]) > na_threshold
  if (all_NAs || t_NA_frac || c_NA_frac) {
    ans <- NA
  } else {
    roc <- ROC::rocdemo.sca(data = x[!is.na(x)], truth = state[!is.na(x)],
      cutpts = seq(0, 100, 1))
    ans <- ROC::AUC(roc)
  }
  return(ans)
}
