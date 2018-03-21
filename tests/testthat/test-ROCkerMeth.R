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
  auc <- single_AUC(c(tumor_toy_table[1,], control_toy_table[1,]), na_threshold=0,
    state = c(rep(T,ncol(tumor_toy_table)), rep(F,ncol(control_toy_table))))
  expect_is(auc, "numeric")
})

test_that("compute_AUC works", {
  expect_error(compute_AUC(1, 1, 100), "Selected")
  expect_is(compute_AUC(head(tumor_toy_table, 10), head(control_toy_table, 10),
      1), "numeric")
})

context("meth_state_finder") ##################################################

test_that("fix_short_segments", {
  x <- c(3,3,3,3,3,3,2,2,2,1,1,2,3,3,3,2,2,1,1,1,1,1,1)
  y <- c(rep(3,6),rep(2,11),rep(1,6))
  expect_equal(fix_short_segments(x, 5), y)
})

test_that("fix_na_segments", {
              # OK                     # too long            # diff flanking sites
  x <- c(1,1, NA, NA, NA, 1,2,2,2,2,3, NA, NA, NA, NA, 3, 1, NA, NA, NA, 2,2,2)
  y <- c(1,1, 1,  1,  1,  1,2,2,2,2,3, NA, NA, NA, NA, 3, 1, NA, NA, NA, 2,2,2)
  expect_equal(fix_na_segments(x, 3), y)
})

test_that("meth_state_finder outputs only NA,1,2,3 states", {
  idx_chr <- which(reference_toy_table[[1]] == "2")
  x <- auc_toy_vector[idx_chr]
  y <- reference_toy_table[[2]][idx_chr]
  meth_states <- meth_state_finder(x, y, 0, 0)
  expect_true(all(meth_states %in% c(NA,1:3)))
})

context("segmentator functions") ######################################
test_that("segmentator returns correct table", {
  idx_chr <- which(reference_toy_table[[1]] == "2")
  meth_states <- meth_state_finder(auc_toy_vector[idx_chr],
    reference_toy_table[[2]][idx_chr], 3, 3)
  tumor_toy_beta_mean <- apply(tumor_toy_table[idx_chr,], 1, mean, na.rm = T)
  normal_beta_mean <- apply(control_toy_table[idx_chr,], 1, mean, na.rm = T)
  single_chr_segs <- segmentator(meth_states, tumor_toy_beta_mean, normal_beta_mean)
  expect_is(single_chr_segs, "data.frame")
  expect_length(single_chr_segs, 4)
  expect_identical(names(single_chr_segs), c("nseg", "state", "avg_beta_diff",
      "pval"))
})

test_that("whole_genome_segmentator works", {
  dmr_table <- whole_genome_segmentator(tumor_toy_table, control_toy_table, auc_toy_vector,
    reference_toy_table)
  expect_is(dmr_table, "data.frame")
  expect_length(dmr_table, 8)
  expect_identical(names(dmr_table),
    c("chr", "start", "end", "nseg", "state", "avg_beta_diff", "pval", "FDR"))
})

context("output") ##############################################################
test_that("compute_z_score works", {
  dmr_table <- whole_genome_segmentator(tumor_toy_table, control_toy_table,
    auc_toy_vector, reference_toy_table)
  sample_score <- compute_z_score(tumor_toy_table, control_toy_table,
    dmr_table, reference_toy_table)
  expect_is(sample_score, "list")
  expect_length(sample_score, 3)
})

context("long test") ########################################################
test_that("TCGA-THCA works", {
  skip("checking")
  data_folder <- "/projects/databases/data/TCGA/clean/harmonized/THCA/"
  skip_if_not(dir.exists(data_folder))
  x <- readRDS(file.path(data_folder, "THCA_DNAm_TP.rds"))
  y <- readRDS(file.path(data_folder, "THCA_DNAm_NT.rds"))
  auc <- readRDS("/projects/packages/ROCkerMeth/data-raw/pancancer_auc.rds")[,"THCA"]
  load("/projects/packages/PAMES/data/illumina450k_hg19.rda")
  idx <- order(illumina450k_hg19$Chromosome, illumina450k_hg19$Genomic_Coordinate)
  x <- x[idx,]
  y <- y[idx,]
  x_mean <- apply(x, 1, mean, na.rm = T)
  y_mean <- apply(y, 1, mean, na.rm = T)
  mean_diff <- x_mean - y_mean
  auc <- auc[idx]
  illumina450k_hg19 <- illumina450k_hg19[idx,]

  dmr_table <- whole_genome_segmentator(x, y, auc, illumina450k_hg19[3:4])
  sample_score <- compute_z_score(x, y, dmr_table, illumina450k_hg19[3:4])
  expect_is(sample_score, "list")
  expect_is(dmr_table, "data.frame")
  expect_is(sample_score, "list")
})
