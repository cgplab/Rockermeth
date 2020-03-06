context("output") ##############################################################
test_that("write_dmr_output works", {
    dmr_table <- suppressMessages(find_dmrs(tumor_toy_table, control_toy_table, auc_toy_vector, reference_toy_table))
    write_dmr(dmr_table, "test")
    expect_true(file.exists("test_hyper.bed"))
    expect_true(file.exists("test_hypo.bed"))
    expect_true(file.exists("test.seg"))
    file.remove("test_hyper.bed", "test_hypo.bed", "test.seg")
})

test_that("write_z_scores works", {
    dmr_table <- suppressMessages(find_dmrs(tumor_toy_table, control_toy_table, auc_toy_vector, reference_toy_table))
    sample_score <- suppressMessages(compute_z_scores(tumor_toy_table, control_toy_table, dmr_table, reference_toy_table))
    write_z_scores(sample_score, "test")
    expect_true(file.exists("test_z_scores.seg"))
    expect_true(file.exists("test_tumor_dmr.seg"))
    expect_true(file.exists("test_control_dmr.seg"))
    file.remove("test_z_scores.seg", "test_tumor_dmr.seg", "test_control_dmr.seg")
})
