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
#' @param max_dist Maximum distance between two contiguous sites in a DMR.
#' @param min_sites Minimum of number of sites to consider a DMR valid.
#'
#' @return A data.frame with information about DMRs
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @keywords internal
segmentator <- function(tumor_beta_mean, control_beta_mean, meth_states, coordinates, max_dist, min_sites){
    assertthat::assert_that(length(tumor_beta_mean) == length(control_beta_mean))
    assertthat::assert_that(length(tumor_beta_mean) == length(meth_states))
    assertthat::assert_that(length(tumor_beta_mean) == length(coordinates))

    beta_diff <- tumor_beta_mean - control_beta_mean

    dmrs <- c()
    # initialize start position and number of sites in dmr
    region_start <- coordinates[1]
    nsites <- 1
    pb <- txtProgressBar(min = 0, max = length(meth_states), style = 3, width=80)
    for (i in seq_along(meth_states)[-1]){ # skip first position
        curr_state <- meth_states[i]
        curr_position <- coordinates[i]
        prev_state <- meth_states[i-1]
        prev_position <- coordinates[i-1]

        # different state or distance threshold define a new region unless in non-diff regions
        if ((curr_state != prev_state || curr_position - prev_position > max_dist && curr_state != 2)){
            region_end <- prev_position
            new_dmr <-
                dplyr::tibble(start = region_start,
                              end   = region_end,
                              nsites,
                              state = prev_state)
            dmrs <- rbind(dmrs, new_dmr)
            # reinitialize variables
            region_start <- curr_position
            nsites <- 1

        } else { # still in the same region
            # increment number of sites
            nsites <- nsites + 1
            # shift end position
        }
        if (i %% 100 == 0)
            setTxtProgressBar(pb, i)
    }
    setTxtProgressBar(pb, i)
    close(pb)
    # add last region
    region_end <- curr_position
    new_dmr <- dplyr::tibble(start = region_start,
                             end   = region_end,
                             nsites = nsites,
                             state = prev_state)
    dmrs <- rbind(dmrs, new_dmr)

    values <- matrix(NA, nrow = nrow(dmrs), ncol = 2,
                     dimnames = list(NULL, c("avg_beta_diff", "p_value")))
    dmrs_idx <- with(dmrs, which(nsites >= min_sites & state != 2))
    start_idx <- c(1, cumsum(dmrs$nsites)[-nrow(dmrs)]+1)
    end_idx   <- cumsum(dmrs$nsites)
    message("Total regions: ", nrow(dmrs))
    message("Total DMRs: ", sum(dmrs$state != 2))
    message("Valid DMRs: ", length(dmrs_idx))
    if (length(dmrs_idx) > 0) {
        pb <- txtProgressBar(min = 0, max = length(dmrs_idx), style = 3, width = 80)
        for (p in seq_along(dmrs_idx)) {
            i <- dmrs_idx[p]
            values[i, 1] <- with(dmrs, mean(beta_diff[start_idx[i]:end_idx[i]], na.rm = TRUE))
            hypothesis <- ifelse(dmrs$state[i] == 3, "greater", "less")
            values[i, 2] <- suppressWarnings(tryCatch(error = function(e) return(NA),
                with(dmrs, wilcox.test(tumor_beta_mean[start_idx[i]:end_idx[i]],
                                       control_beta_mean[start_idx[i]:end_idx[i]],
                                       hypothesis)[["p.value"]])))
            setTxtProgressBar(pb, p)
        }
        close(pb)
    }
    dmrs <- dplyr::bind_cols(dmrs, dplyr::as_tibble(values))
    return(dmrs)
}
