context("segmentator functions") ######################################
test_that("segmentator returns correct table", {
    idx_chr <- with(reference_toy_table, which(Chromosome == "1" & Chromosome_Arm == "p"))
    auc_sd <- sd(auc_toy_vector, na.rm = TRUE)
    auc <- auc_toy_vector[idx_chr]
    coordinates <- reference_toy_table$Genomic_Coordinate[idx_chr]
    idx_not_NA <- which(!is.na(auc))
    meth_states <- meth_state_finder(auc[idx_not_NA], coordinates[idx_not_NA], auc_sd,
                                     pt_start = 0.05, normdist = 1e5, ratiosd = 0.4, mu=.1, use_trunc=FALSE)

    tumor_toy_beta_mean <- apply(tumor_toy_table[idx_chr,], 1, mean, na.rm = TRUE)
    control_beta_mean <- apply(control_toy_table[idx_chr,], 1, mean, na.rm = TRUE)
    dmrs <- segmentator(tumor_toy_beta_mean[idx_not_NA],
                        control_beta_mean[idx_not_NA],
                        meth_states,
                        coordinates[idx_not_NA], Inf)

    expect_is(dmrs, "data.frame")
    expect_identical(names(dmrs), c("start", "end", "nsites", "state", "mean_beta_diff", "mean_beta_tumor", "mean_beta_control", "p_value"))
})

test_that("find_dmrs works", {
    dmr_table <- find_dmrs(tumor_toy_table, control_toy_table, auc_toy_vector, reference_toy_table)
    # head(as.data.frame(dmr_table))
    #   chr  start    end nsites state mean_beta_diff mean_beta_tumor mean_beta_control      p_value      fdr
    # 1   1  15865 763127     47     2     -0.1250561        28.09480          28.21986           NA           NA
    # 2   1 779995 876551     94     1    -14.4714773        36.97280          51.44428 0.0000664017 0.0006451733
    # 3   1 876870 901062     51     2     -0.8461393        56.63192          57.47806           NA           NA
    # 4   1 901449 903543     10     1     -8.4722725        10.96582          19.43810 0.0715700708 0.1283263863
    # 5   1 905988 910223      7     2      2.4901404        59.74864          57.25850           NA           NA
    # 6   1 911995 937048     38     1    -10.6365439        21.95869          32.59524 0.0627186188 0.1217088467
    expect_is(dmr_table, "data.frame")
    expect_length(dmr_table, 10)
    expect_identical(names(dmr_table), c("chr", "start", "end", "nsites", "state", "mean_beta_diff", "mean_beta_tumor", "mean_beta_control", "p_value", "fdr"))
})
