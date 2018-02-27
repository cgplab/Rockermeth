#' Compute Area Under Curve for a matrix of samples
#'
#' This function compute the Area Under Curve that is then used to define the segregation
#' between tumor and normal samples accordingly to their methylation (beta) values.
#'
#' @param tumor_table A matrix of methylation (beta) values produced by a 450k Illumina
#' BeadChip from a \bold{tumor} sample
#' @param control_table A matrix of methylation (beta) values produced by a 450k Illumina
#' BeadChip from a \bold{control} sample
#' @param nclust Number of clusters to use in parallel
#' @param NA_thr Fraction of NAs (considered independently in tumor and
#' control samples) above which a site will not be selected (default=0)
#' @return A vector of AUC scores
#' @export
compute_AUC <- function(tumor_table, control_table, nclust = 1, NA_thr = 0) {
  nclust <- as.integer(nclust)
  if (nclust > parallel::detectCores()){
    stop(sprintf("Selected %i cores but system has %i cores.", nclust, parallel::detectCores()))
  }

  beta_table <- as.matrix(cbind(tumor_table, control_table))
  sample_state <- c(rep(T, ncol(tumor_table)), rep(F, ncol(control_table)))
  if (all(beta_table >= 0, na.rm = T) && all(beta_table <= 1, na.rm = T) && is.double(beta_table[1,1])) {
    beta_table <- beta_table*100
    storage.mode(beta_table) <- "integer"
  } else {
    stop("tumor_table and control_table must have fraction values")
  }

  NA_thr <- as.numeric(NA_thr)
  stopifnot(NA_thr >= 0 || NA_thr < 1)

  cl <- parallel::makeCluster(nclust)
  auc <- parallel::parApply(cl, beta_table, 1, single_AUC, 
    state = sample_state, NA_thr = NA_thr)
  parallel::stopCluster(cl)
  return(auc)
}

#' Compute AUC for a single vector
#'
#' Compute AUC returning NA if NA samples are more than threshold
#' @param x integer vector (range 1-100)
#' @param state logical vector
#' @param NA_thr numeric value: if fraction of NAs is higher than threshold,
#' either in tumor or control samples, return NA
#' @keywords internal
single_AUC <- function(x, state, NA_thr) {
  all_NAs <- all(is.na(x)) 
  t_NA_frac <- sum(is.na(x[state])) / length(x[state]) > NA_thr 
  c_NA_frac <- sum(is.na(x[!state])) / length(x[!state]) > NA_thr 
  if (all_NAs || t_NA_frac || c_NA_frac) {
    ans <- NA
  } else {
    roc <- ROC::rocdemo.sca(data = x[!is.na(x)], truth = state[!is.na(x)],
      cutpts = seq(0, 100, 1))
    ans <- ROC::AUC(roc)
  }
  return(ans)
}
