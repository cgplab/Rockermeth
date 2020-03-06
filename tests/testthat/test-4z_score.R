context("z-score") ##############################################################
test_that("compute_z_scores 'default' works", {
    dmr_table <- find_dmrs(tumor_toy_table, control_toy_table, auc_toy_vector, reference_toy_table)
    sample_score <- compute_z_scores(tumor_toy_table, control_toy_table, dmr_table, reference_toy_table)
    expect_is(sample_score, "list")
    expect_length(sample_score, 4)
})

test_that("compute_z_scores 'custom' works", {
    dmr_table <- find_dmrs(tumor_toy_table, control_toy_table, auc_toy_vector, reference_toy_table)
    dmr_table_r <- dplyr::filter(dmr_table, q_value < 0.05, state == 3) %>% head(20)
    sample_score_r <- compute_z_scores(tumor_toy_table, control_toy_table, dmr_table_r, reference_toy_table, method = "custom", min_sites = 10)
    expect_is(sample_score_r, "list")
    expect_length(sample_score_r, 4)
    expect_equal(dim(sample_score_r[[1]]), c(20,79))
})
