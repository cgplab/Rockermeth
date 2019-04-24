#' Define Differentially Methylated Regions
#'
#' Given a set of methylation states produced by [meth_states_finder] and a set
#' of methylation (beta) values differences between tumor and control samples,
#' this function defines stretches of equal states (segments).
#'
#' @param tumor_beta_mean A vector of average methylation (beta) values.
#' @param control_beta_mean A vector of average methylation (beta) values.
#' @param meth_states A vector of methylation states (1, 2, 3, NA).
#' @param coordinates Ordered genomic positions of methylation sites
#' @param max_dist Maximum distance between sites within same DMRs
#'
#' @return A data.frame with information about DMRs
#' @export
#' @keywords internal
segmentator <- function(tumor_beta_mean, control_beta_mean, meth_states, coordinates, max_dist){
  assertthat::assert_that(length(tumor_beta_mean) == length(control_beta_mean))
  assertthat::assert_that(length(tumor_beta_mean) == length(meth_states))
  assertthat::assert_that(length(tumor_beta_mean) == length(coordinates))

  # insert cuts (0s) between sites more distant than max_dist
  # clever solution from https://stackoverflow.com/a/1495204/3243916
  distant_sites <- which(c(diff(coordinates), 0) > max_dist)
  split_states <- c(meth_states, rep(0, length(distant_sites)))
  idx <- c(seq_along(meth_states), distant_sites + 0.5)
  split_states <- split_states[order(idx)]
  meth_states <- split_states

  beta_diff <- tumor_beta_mean - control_beta_mean

  rle_states <- S4Vectors::Rle(meth_states)
  rle_value  <- S4Vectors::runValue(rle_states)
  idx <- which(S4Vectors::runValue(rle_states) != 0)
  rle_value  <- rle_value[idx]
  rle_length <- S4Vectors::runLength(rle_states)[idx]
  rle_n      <- length(rle_value)
  rle_start  <- S4Vectors::start(rle_states)[idx]
  rle_end    <- S4Vectors::end(rle_states)[idx]

  values <- matrix(NA, nrow = rle_n, ncol = 2, dimnames = list(NULL, c("avg_beta_diff", "p_value")))

  for (i in seq_len(rle_n)) {
    values[i, 1] <- mean(beta_diff[rle_start[i]: rle_end[i]], na.rm = TRUE)

    if (rle_value[i] != 2) {
      hypothesis <- dplyr::if_else(rle_value[i] == 3, "greater", "less")
      values[i, 2] <- suppressWarnings(
        tryCatch(error = function(e) return(NA),
            wilcox.test(tumor_beta_mean[rle_start[i]:rle_end[i]],
                        control_beta_mean[rle_start[i]:rle_end[i]],
                        hypothesis)[["p.value"]]))
    }
  }
  dmrs <- dplyr::tibble(start = coordinates[rle_start],
                        end   = coordinates[rle_end],
                        nsites = rle_length,
                        state = rle_value)
  dmrs <- dplyr::bind_cols(dmrs, dplyr::as_tibble(values))
  return(dmrs)
}
