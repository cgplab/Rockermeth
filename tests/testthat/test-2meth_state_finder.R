context("meth_state_finder") ##################################################
test_that("meth_state_finder works", {
    idx_chr <- with(reference_toy_table, which(Chromosome == "1" & Chromosome_Arm == "p"))
    auc_sd <- sd(auc_toy_vector, na.rm = TRUE)
    auc <- auc_toy_vector[idx_chr]
    coordinates <- reference_toy_table$Genomic_Coordinate[idx_chr]
    idx_not_NA <- which(!is.na(auc))
    meth_states <- meth_state_finder(auc[idx_not_NA], coordinates[idx_not_NA], auc_sd,
                                     pt_start = 0.05, normdist = 1e5, ratiosd = 0.4, mu=.1, use_trunc=FALSE)
    expect_equal(as.numeric(table(meth_states)), c(3329, 18566, 775))
})
