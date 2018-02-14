context("AUC functions") #####################################################

test_that("single_AUC returns correct value", {
  a <- c(runif(10, 51, 100), rep(NA, 2))
  b <- c(runif(10,  0,  50), rep(NA, 2))
  d <- c(runif(10, 51, 100), rep(NA, 1000))
  e <- c(runif(10,  0,  50), rep(NA, 1000))
  auc <- single_AUC(c(a, b), c(rep(T, length(a)), rep(F, length(b))), .2)
  expect_equal(auc, 1)
  auc <- single_AUC(c(a, b), c(rep(T, length(a)), rep(F, length(b))), 0)
  expect_equal(auc, NA)
  auc <- single_AUC(c(d, e), c(rep(T, length(d)), rep(F, length(e))), .9)
  expect_equal(auc, NA)
  auc <- single_AUC(c(tumor_table[1,], control_table[1,]),
    state= c(rep(T,ncol(tumor_table)), rep(F,ncol(control_table))), NA_thr = 0)
  expect_is(auc, "numeric")
})

test_that("compute_AUC works", {
  expect_error(compute_AUC(1, 1, 100), "Selected")
  expect_is(compute_AUC(head(tumor_table, 10), head(control_table, 10), 1), "numeric")
})

context("meth_state_finder") ##################################################

test_that("fix_short_segments corrects short methylated streches", {
  x <- c(3,3,3,3,3,3,2,2,2,1,1,2,3,3,3,2,2,1,1,1,1,1,1)
  expect_equal(fix_short_segments(x, 5), c(rep(3,6),rep(2,11),rep(1,6)))
})

test_that("fix_na_states set NA stretches to 2", {
  x <- c(2,2,3,3, NA, NA, NA, 2,2,2)
  y <- c(2,2,3,3, 3,  3,  3, 2,2,2)
  expect_equal(fix_na_states(x, 4), y)

  x <- c(2,2,3,3, NA, NA, NA, 2,2,2)
  y <- c(2,2,3,3, NA, NA, NA, 2,2,2)
  expect_equal(fix_na_states(x, 3), y)

  x <- c(NA, NA, NA, 1,1,1,1, 3,3,3)
  y <- c(1,  1,  1,  1,1,1,1, 3,3,3)
  expect_equal(fix_na_states(x, 4), y)

  x <- c(NA, NA, NA, 1,1,1,1, 3,3,3)
  y <- c(NA, NA, NA, 1,1,1,1, 3,3,3)
  expect_equal(fix_na_states(x, 3), y)

  x <- c(2,2,2, 1,1,1,1, NA, NA, NA)
  y <- c(2,2,2, 1,1,1,1, 1,  1,  1)
  expect_equal(fix_na_states(x, 4), y)

  x <- c(2,2,2, 1,1,1,1, NA, NA, NA)
  y <- c(2,2,2, 1,1,1,1, NA, NA, NA)
  expect_equal(fix_na_states(x, 3), y)

  x <- c(1,1, NA, NA, NA, NA, 1,2,2,2,2,3, NA, NA, NA, 1, 1, NA, NA, NA, NA, 2,2,2)
  y <- c(1,1, 1,  1,  1,  1,  1,2,2,2,2,3, 3,  3,  3,  1, 1, NA, NA, NA, NA, 2,2,2)
  expect_equal(fix_na_states(x, 4), y)
})

test_that("meth_state_finder outputs only NA,1,2,3 states", {
  idx_chr <- which(reference_table[[1]] == "2")
  meth_states <- meth_state_finder(auc_vector[idx_chr], reference_table[[2]][idx_chr])
  expect_true(all(meth_states %in% c(NA, 1:3)))
})

context("segmentator functions") ######################################
test_that("segmentator returns correct table", {
  idx_chr <- which(reference_table[[1]] == "2")
  meth_states <- meth_state_finder(auc_vector[idx_chr], reference_table[[2]][idx_chr])
  tumor_beta_mean <- apply(tumor_table[idx_chr,], 1, mean, na.rm = T)
  normal_beta_mean <- apply(control_table[idx_chr,], 1, mean, na.rm = T)
  single_chr_segs <- segmentator(meth_states, tumor_beta_mean, normal_beta_mean)
  expect_is(single_chr_segs, "data.frame")
  expect_length(single_chr_segs, 4)
  expect_identical(names(single_chr_segs), c("nseg", "state", "avg_beta_diff", "pval"))
})

test_that("whole_genome_segmentator works", {
  all_breaks <- whole_genome_segmentator(tumor_table, control_table, auc_vector,
    reference_table)
  expect_is(all_breaks, "data.frame")
  expect_length(all_breaks, 8)
  expect_identical(names(all_breaks),
    c("chr", "start", "end", "nseg", "state", "avg_beta_diff", "pval", "FDR"))
})

context("output") ##############################################################
test_that("calc_sample_score works", {
  dmr_table <- whole_genome_segmentator(tumor_table, control_table, auc_vector, reference_table)
  sample_score <- calc_sample_score(tumor_table, control_table, dmr_table, reference_table)
  expect_is(sample_score, "list")
  expect_length(sample_score, 3)
})
