#' Compute a Z-like score to estimate quality of DMRs
#'
#' Assayed CpG may differ from the sites used in auc computation and
#' segmentation. However, using the same assay for 'Normal' and 'Tumor' is
#' strongly recommended.
#'
#' @param tumor_table A matrix of beta-values (fraction) from tumor samples.
#' @param control_table A matrix of beta-values (fraction) from normal/control samples.
#' @param dmr_table A data.frame reporting the genomic location of DMRs
#' (chromosome, start, end) (likely produced by [find_dmrs]).
#' @param reference_table A data.frame reporting the genomic coordinates of
#' each CpG site in tumor and control matrices.
#' @param min_size Minimum number of CpG sites inside DMR to compute Z-score
#' (default = 3, return NA if lower).
#' @return A list of 4 tables: z-scores of DMRs, median beta of DMRs in
#' tumor samples, median beta of DMRs in normal/control samples and fraction of
#' NA CpG sites within DMRs.
#'
#' @importFrom stats mad median
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
compute_z_scores <- function(tumor_table, control_table, dmr_table, reference_table, min_size=3) {

    # check parameters
    assertthat::assert_that(is.data.frame(reference_table))
    assertthat::assert_that(length(reference_table) >= 2)
    assertthat::assert_that(nrow(tumor_table) == nrow(control_table))
    assertthat::assert_that(nrow(tumor_table) == nrow(reference_table))
    assertthat::assert_that(is.data.frame(dmr_table))
    assertthat::assert_that(length(dmr_table) >= 3)

    # check chromosome names
    assertthat::assert_that(length(intersect(reference_table[[1]], dmr_table[[1]])) > 0,
                            msg="No shared chromosomes between 'reference_table' and 'dmr_table.' Check chromosome names.")

    beta_table <- as.matrix(cbind(tumor_table, control_table))
    diff_range <- diff(range(beta_table, na.rm = TRUE))
    assertthat::assert_that(diff_range > 1, diff_range <= 100, msg = "For computation efficiency, convert tumor table to percentage values.")
    beta_table <- round(beta_table)
    storage.mode(beta_table) <- "integer"

    # compute z-scores
    sample_state <- c(rep(TRUE, ncol(tumor_table)), rep(FALSE, ncol(control_table)))
    tumor_dmr_beta   <- matrix(NA, nrow(dmr_table), ncol(tumor_table))
    control_dmr_beta <- matrix(NA, nrow(dmr_table), ncol(control_table))
    z_scores         <- matrix(NA, nrow(dmr_table), ncol(tumor_table))
    na_frac          <- matrix(NA, nrow(dmr_table), ncol(tumor_table)) # NA fraction matrix as a quality feedback

    # use CpG sites located within DMRs
    sites <- GenomicRanges::GRanges(seqnames = reference_table[[1]],
                                    ranges = IRanges::IRanges(start = reference_table[[2]], width = 1),
                                    idx = seq_len(nrow(reference_table)))
    dmrs <- GenomicRanges::GRanges(seqnames = dmr_table[[1]],
                                   ranges = IRanges::IRanges(start = dmr_table[[2]], end = dmr_table[[3]]))
    overlaps <- as.data.frame(GenomicRanges::findOverlaps(sites, dmrs))
    dmr_idxs <- unique(overlaps$subjectHits)

    insuff_segs <- 0
    message(sprintf("[%s] Computing DMR median beta", Sys.time()))
    pb <- txtProgressBar(min = 0, max = length(dmr_idxs), style = 3, width=80)
    for (i in seq_along(dmr_idxs)) {
        idx <- dmr_idxs[i]
        if (sum(overlaps$subjectHits == idx) >= min_size) {
            idx_dmr <- overlaps$queryHits[overlaps$subjectHits == idx]
            tumor_dmr_beta[idx,] <-
                apply(beta_table[idx_dmr, sample_state, drop = FALSE], 2, median, na.rm = TRUE)
            control_dmr_beta[idx,] <-
                apply(beta_table[idx_dmr, !sample_state, drop = FALSE], 2, median, na.rm = TRUE)
            ## Compute percentage of NA values for each DMR, to get a feedback on how reliable the result is
            na_frac[idx,] <-
                apply(beta_table[idx_dmr, sample_state, drop = FALSE], 2, function(x) sum(is.na(x))/length(x))
        } else {
            insuff_segs <- insuff_segs + 1
        }
        if (i %% 100 == 0){
            setTxtProgressBar(pb, i)
        }
    }
    setTxtProgressBar(pb, i)
    close(pb)

    dmr_without_signal <- nrow(dmr_table) - length(dmr_idxs)
    message(sprintf("DMRs with no probes: %i", dmr_without_signal))
    message(sprintf("DMRs with not enough probes: %i ", insuff_segs))
    message(sprintf("[%s] Computing z-scores", Sys.time()))

    control_median <- apply(control_dmr_beta, 1, median, na.rm=TRUE)
    control_median_abs_dev <- apply(control_dmr_beta, 1, mad, na.rm=TRUE)
    control_median_abs_dev[dplyr::between(control_median_abs_dev, 0, 1)] <- 1
    z_scores <- (tumor_dmr_beta - control_median) / control_median_abs_dev

    rnames <- sprintf("chr%s:%s-%s", dmr_table[[1]], dmr_table[[2]], dmr_table[[3]])

    dimnames(tumor_dmr_beta)   <- list(rnames, colnames(beta_table)[sample_state])
    dimnames(control_dmr_beta) <- list(rnames, colnames(beta_table)[!sample_state])
    dimnames(z_scores)         <- list(rnames, colnames(beta_table)[sample_state])
    dimnames(na_frac)          <- list(rnames, colnames(beta_table)[sample_state])

    return(list(z_scores = z_scores,
                tumor_dmr_beta = tumor_dmr_beta,
                control_dmr_beta = control_dmr_beta,
                na_frac = na_frac))
}
