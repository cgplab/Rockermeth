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
  # check parameters
  ncores <- as.integer(ncores)
  system_cores <- parallel::detectCores()
  assertthat::assert_that(ncores < system_cores)

  na_threshold <- as.numeric(na_threshold)
  assertthat::assert_that(na_threshold >= 0 || na_threshold < 1)

  beta_table <- as.matrix(cbind(tumor_table, control_table))
  diff_range <- diff(range(beta_table, na.rm = TRUE))
  if (diff_range <= 1 || diff_range > 100) {
    stop(paste0("For computation efficiency please convert tumor and control",
        "tables to percentage value."))
  } else {
    beta_table <- round(beta_table)
    storage.mode(beta_table) <- "integer"
  }
  sample_state <- c(rep(TRUE, ncol(tumor_table)), rep(FALSE, ncol(control_table)))

  cl <- parallel::makeCluster(ncores)
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
  assertthat::assert_that(is.integer(x))
  all_NA <- all(is.na(x))
  tumor_NA <- sum(is.na(x[state])) / length(x[state]) > na_threshold
  control_NA <- sum(is.na(x[!state])) / length(x[!state]) > na_threshold
  if (all_NA || tumor_NA || control_NA ) {
    return(NA)
  } else {
    roc <- ROC::rocdemo.sca(state[!is.na(x)], x[!is.na(x)],
      cutpts = seq(0, 100, 1))
    return(ROC::AUC(roc))
  }
}
