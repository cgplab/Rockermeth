#' Compute a Z-like score to estimate quality of DMRs
#'
#' Assayed CpG may differ from the sites used in AUC computation and
#' segmentation. However, using the same assay for 'Normal' and 'Tumor' is
#' strongly recommended.
#'
#' @param tumor_table A matrix of beta-values (fraction) from tumor samples.
#' @param control_table A matrix of beta-values (fraction) from normal/control samples.
#' @param dmr_table A data.frame reporting the genomic location of DMRs
#' (chromosome, start, end) (likely produced by [find_dmrs]).
#' @param reference_table A data.frame reporting the genomic coordinates of
#' each CpG site in tumor and control matrices.
#' @param method "default" reports only statistically significant DMRs: used
#' for a DMR set generated with [find_dmrs]; "custom" compute Z-scores of regions
#' covered by a minimum number of CpG sites: used to compare regions obtained with
#' different tools
#' @param q_value_thr Minimum q-value of a DMR to
#' compute a Z-score (used only in "default" analysis; default = 0.05).
#' @param min_sites Minimum required number of CpG sites within a DMR to
#' compute a Z-score (used only in "custom" analysis; default = 5).
#' @param ncores Number of parallel processes to use for parallel computing.
#' @return A list of 4 tables: z-scores of DMRs, median beta of DMRs in
#' tumor samples, median beta of DMRs in normal/control samples and fraction of
#' NA CpG sites within DMRs.
#' @examples
#' auc <- compute_AUC(tumor_example, control_example)
#' dmr_set <- find_dmrs(tumor_example, control_example, auc, reference_example, min_sites = 10)
#' compute_z_scores(tumor_example, control_example, dmr_set, reference_example)
#' @importFrom stats mad median
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
compute_z_scores <- function(tumor_table, control_table, dmr_table,
                             reference_table, method=c("default", "custom"),
                             q_value_thr = 0.05, min_sites=5, ncores=1) {
    message(sprintf("[%s] Z-scores analysis", Sys.time()))
    # check parameters
    system_cores <- parallel::detectCores()
    
    if (!is.na(system_cores)){
        assertthat::assert_that(ncores < system_cores)
        }
    
    method <- match.arg(method)
    assertthat::assert_that(is.data.frame(reference_table))
    assertthat::assert_that(length(reference_table) >= 2)
    assertthat::assert_that(nrow(tumor_table) == nrow(control_table))
    assertthat::assert_that(nrow(tumor_table) == nrow(reference_table))
    # check chromosome names
    assertthat::assert_that(length(intersect(reference_table[[1]], dmr_table[[1]])) > 0,
                            msg="No shared chromosomes between 'reference_table' and 'dmr_table.' Check chromosome's names.")

    beta_table <- as.matrix(cbind(tumor_table, control_table))
    beta_table <- round(beta_table)
    storage.mode(beta_table) <- "integer"

    if (method == "default"){
        dmr_table <- dplyr::filter(dmr_table, q_value < q_value_thr)
        if (nrow(dmr_table) == 0){
            warning("No DMRs with significant q-value retrieved: no output produced.")
            return(NULL)
        }
    }

    tumor_dmr_beta   <- matrix(NA, nrow(dmr_table), ncol(tumor_table))
    control_dmr_beta <- matrix(NA, nrow(dmr_table), ncol(control_table))
    z_scores         <- matrix(NA, nrow(dmr_table), ncol(tumor_table))
    na_frac          <- matrix(NA, nrow(dmr_table), ncol(beta_table)) # NA fraction matrix as a quality feedback

    rnames <- sprintf("%s:%i-%i", dmr_table[[1]], dmr_table[[2]], dmr_table[[3]])
    dimnames(tumor_dmr_beta)   <- list(rnames, colnames(tumor_table))
    dimnames(control_dmr_beta) <- list(rnames, colnames(control_table))
    dimnames(z_scores)         <- list(rnames, colnames(tumor_table))
    dimnames(na_frac)          <- list(rnames, colnames(beta_table))

    # find CpG sites located within DMRs
    message(sprintf("[%s] Find overlaps", Sys.time()))
    dmrs_ranges <- GenomicRanges::GRanges(seqnames = dmr_table[[1]],
                                          ranges = IRanges::IRanges(start = dmr_table[[2]],
                                                                    end = dmr_table[[3]]))
    sites_ranges <- GenomicRanges::GRanges(seqnames = reference_table[[1]],
                                           ranges = IRanges::IRanges(start = reference_table[[2]],
                                                                     width = 1))
    overlaps <- GenomicRanges::findOverlaps(dmrs_ranges, sites_ranges)
    # split sites indexes by DMRs
    sites_idx_list <- split(S4Vectors::subjectHits(overlaps),
                            S4Vectors::queryHits(overlaps))
    if (method == "custom"){
        dmrs_nsites <- as.vector(sapply(sites_idx_list, length))
        valid_dmrs <- which(dmrs_nsites >= min_sites)
    } else {
        valid_dmrs <- seq_along(sites_idx_list)
    }

    message(sprintf("[%s] Compute statistics", Sys.time()))
    ## Compute median beta and percentage of NA values per DMR per sample
    dmrs_info <- parallel::mclapply(mc.cores = ncores, sites_idx_list[valid_dmrs], function(idx) {
        y <- apply(beta_table[idx,, drop = FALSE], 2, function(x) {
                dmr_beta <- median(x, na.rm = TRUE)
                na_frac <- sum(is.na(x))/length(x)
                return(cbind(dmr_beta, na_frac))
        })
        return(y)
    })

    dmrs_idx <- as.integer(names(sites_idx_list[valid_dmrs]))
    sample_state <- c(rep(TRUE, ncol(tumor_table)), rep(FALSE, ncol(control_table)))
    tumor_dmr_beta[dmrs_idx,]     <- t(sapply(dmrs_info, function(x) x[1, sample_state]))
    control_dmr_beta[dmrs_idx,]   <- t(sapply(dmrs_info, function(x) x[1, !sample_state]))
    na_frac[dmrs_idx,]            <- t(sapply(dmrs_info, function(x) x[2,]))

    cl <- parallel::makeCluster(ncores)
    control_samples_statistics <- t(parallel::parApply(cl, control_dmr_beta, 1, function(x) {
        y <- c(Median=median(x, na.rm=TRUE), MAD=mad(x, na.rm=TRUE))
        y[2] <- dplyr::if_else(dplyr::between(y[2], 0, 1), 1, y[2]) # set low MAD to 1
        return(y)
    }))
    parallel::stopCluster(cl)
    z_scores <- sweep(tumor_dmr_beta, 1, control_samples_statistics[,1])
    z_scores <- sweep(z_scores, 1, control_samples_statistics[,2], "/")
    if (method == "custom") {
        message(sprintf("* DMRs with insufficient number of sites: %i/%i",
                        sum(dmrs_nsites < min_sites), nrow(dmr_table)))
    }

    message(sprintf("[%s] Done", Sys.time()))
    return(list(z_scores = z_scores,
                tumor_dmr_beta = tumor_dmr_beta,
                control_dmr_beta = control_dmr_beta,
                NA_frac = na_frac))
}
