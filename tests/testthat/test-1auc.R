context("AUC functions") #####################################################

test_that("single_AUC returns correct values", {
    a <- as.integer(c(runif(10, 51, 100), rep(NA, 2)))
    b <- as.integer(c(runif(10,  0,  50), rep(NA, 2)))
    d <- as.integer(c(runif(10, 51, 100), rep(NA, 1000)))
    e <- as.integer(c(runif(10,  0,  50), rep(NA, 1000)))
    auc <- single_AUC(c(a, b), c(rep(TRUE, length(a)), rep(FALSE, length(b))))
    expect_equal(auc, 1)
    auc <- single_AUC(c(tumor_toy_table[1,], control_toy_table[1,]),
                      c(rep(TRUE, ncol(tumor_toy_table)), rep(FALSE, ncol(control_toy_table))))
    expect_is(auc, "numeric")
})

test_that("compute_AUC works", {
    expect_error(compute_AUC(tumor_toy_table/100, control_toy_table/100), "For computation efficiency")
    expect_warning(compute_AUC(head(tumor_toy_table, 10), head(control_toy_table, 10), na_threshold = .4), "deprecated")
    auc <- compute_AUC(head(tumor_toy_table, 10), head(control_toy_table, 10))
    expect_is(auc, "numeric")
    auc <- compute_AUC(head(tumor_toy_table, 10), head(control_toy_table, 10), simplify = FALSE)
    expect_is(auc, "data.frame")
})
