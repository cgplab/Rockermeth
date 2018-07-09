#' Define Differentially Methylated Regions
#'
#' Given a set of methylation states produced by [meth_states_finder] and a set of
#' methylation (beta) values differences between tumor and control samples, this fucntion
#' defines stretches of equal states (segments).
#'
#' @param meth_states a vector of methylation states (1, 2, 3, NA)
#' @param tumor_beta_mean a vector of average methylation (beta) values
#' @param control_beta_mean a vector of average methylation (beta) values
#'
#' @return a data.frame with information about segments
#' @keywords internal
segmentator <- function(meth_states, tumor_beta_mean, control_beta_mean){
  assertthat::assert_that(length(tumor_beta_mean) == length(control_beta_mean))
  beta_diff <- tumor_beta_mean - control_beta_mean

  rle_states <- S4Vectors::Rle(meth_states)
  rle_length <- S4Vectors::runLength(rle_states)
  rle_value  <- S4Vectors::runValue(rle_states)
  rle_n      <- S4Vectors::nrun(rle_states)
  rle_start  <- S4Vectors::start(rle_states)
  rle_end    <- S4Vectors::end(rle_states)

  values <- matrix(NA, rle_n, 2,
    dimnames = list(NULL, c("avg_beta_diff", "p_value")))

  for (i in seq_len(rle_n)) {
    values[i, 1] <- mean(beta_diff[rle_start[i]: rle_end[i]], na.rm = TRUE)

    if (rle_value[i] != 2) {
      hypothesis <- dplyr::if_else(rle_value[i] == 3, "greater", "less")
      values[i, 2] <- suppressWarnings(tryCatch(wilcox.test(
        tumor_beta_mean[rle_start[i]:rle_end[i]],
        control_beta_mean[rle_start[i]:rle_end[i]], hypothesis)$p.value,
        error = function(e) NA))
    }
  }
  return(data.frame(nseg = rle_length, state = rle_value, values))
}
