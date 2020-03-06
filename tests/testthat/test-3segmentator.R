cjntext("segmentator functions") ######################################
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
    expect_identical(names(dmrs), c("start", "end", "nsites", "state", "mean_beta_diff", "p_value"))
})

test_that("find_dmrs works", {
    dmr_table <- find_dmrs(tumor_toy_table, control_toy_table, auc_toy_vector, reference_toy_table)
    # head(dmr_table)
    #   chr   start     end nsites state mean_beta_diff   p_value   q_value
    # 1   1   15865  991567    351     2     -2.706297        NA        NA
    # 2   1  994498  997269     10     1     -6.432670 0.2179361 0.3007743
    # 3   1  997858 1053027     91     2     -2.893905        NA        NA
    # 4   1 1053794 1067099     20     1     -7.615732 0.1896548 0.3007743
    # 5   1 1067223 1079622     23     2     -3.314123        NA        NA
    # 6   1 1079879 1093335     10     1     -9.431163 0.2406255 0.3007743
    expect_is(dmr_table, "data.frame")
    expect_length(dmr_table, 8)
    expect_identical(names(dmr_table),
                     c("chr", "start", "end", "nsites", "state", "mean_beta_diff", "p_value", "q_value"))
})
