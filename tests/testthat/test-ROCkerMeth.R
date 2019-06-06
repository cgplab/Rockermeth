context("AUC functions") #####################################################
a <- as.integer(c(runif(10, 51, 100), NA))
b <- as.integer(c(runif(10,  0,  50), NA))
s <- c(rep(TRUE, length(a)), rep(FALSE, length(b)))
single_AUC(c(a, b), s)

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
  expect_error(compute_AUC(1, 1, 100), "ncores not less")
  expect_error(compute_AUC(tumor_toy_table/100, control_toy_table/100),
    "For computation efficiency")
  expect_is(compute_AUC(head(tumor_toy_table, 10), head(control_toy_table, 10)),
      "numeric")
})

context("meth_state_finder") ##################################################

test_that("meth_state_finder works", {
  idx_chr <- which(reference_toy_table[[1]] == "2")
  auc_sd <- sd(auc_toy_vector, na.rm = TRUE)
  auc <- auc_toy_vector[idx_chr]
  coordinates <- reference_toy_table[[2]][idx_chr]
  idx_not_NA <- which(!is.na(auc))
  meth_states <- meth_state_finder(auc[idx_not_NA], coordinates[idx_not_NA], auc_sd,
    pt_start = 0.05, normdist = 1e5, ratiosd = 0.4, mu=.1, use_trunc=FALSE)
  # table(meth_states)
  expect_equal(as.numeric(table(meth_states)), c(743, 27170, 31))
})

context("segmentator functions") ######################################
test_that("segmentator returns correct table", {
  idx_chr <- which(reference_toy_table[[1]] == "2")
  auc_sd <- sd(auc_toy_vector, na.rm = TRUE)
  auc <- auc_toy_vector[idx_chr]
  coordinates <- reference_toy_table[[2]][idx_chr]
  idx_not_NA <- which(!is.na(auc))
  meth_states <- meth_state_finder(auc[idx_not_NA], coordinates[idx_not_NA], auc_sd,
    pt_start = 0.05, normdist = 1e5, ratiosd = 0.4, mu=.25, use_trunc=FALSE)

  tumor_toy_beta_mean <- apply(tumor_toy_table[idx_chr,], 1, mean, na.rm = TRUE)
  control_beta_mean <- apply(control_toy_table[idx_chr,], 1, mean, na.rm = TRUE)
  dmrs <- segmentator(tumor_toy_beta_mean[idx_not_NA], control_beta_mean[idx_not_NA], meth_states, coordinates[idx_not_NA], Inf, 3)

  expect_is(dmrs, "data.frame")
  expect_identical(names(dmrs), c("start", "end", "nsites", "state", "avg_beta_diff", "p_value"))
})

test_that("find_dmrs works", {
  dmr_table <- find_dmrs(tumor_toy_table, control_toy_table,
                                        auc_toy_vector, reference_toy_table)
  # head(dmr_table)
  #   chr   start     end nsites state avg_beta_diff   p_value   q_value
  # 1   1   15865  991567    351     2     -2.706297        NA        NA
  # 2   1  994498  997269     10     1     -6.432670 0.2179361 0.3007743
  # 3   1  997858 1053027     91     2     -2.893905        NA        NA
  # 4   1 1053794 1067099     20     1     -7.615732 0.1896548 0.3007743
  # 5   1 1067223 1079622     23     2     -3.314123        NA        NA
  # 6   1 1079879 1093335     10     1     -9.431163 0.2406255 0.3007743
  expect_is(dmr_table, "data.frame")
  expect_length(dmr_table, 8)
  expect_identical(names(dmr_table),
    c("chr", "start", "end", "nsites", "state", "avg_beta_diff", "p_value", "q_value"))
})

context("ouput") ##############################################################
test_that("compute_z_scores and write_output works", {
  dmr_table <- find_dmrs(tumor_toy_table, control_toy_table,
    auc_toy_vector, reference_toy_table)
  sample_score <- compute_z_scores(tumor_toy_table, control_toy_table,
    dmr_table, reference_toy_table, 1)
  expect_is(sample_score, "list")
  expect_length(sample_score, 4)
  write_output("test", dmr_table, sample_score, 0.8)
  expect_true(file.exists("test_hyper.bed"))
  expect_true(file.exists("test_hypo.bed"))
  expect_true(file.exists("test.seg"))
  expect_true(file.exists("test_z_scores.seg"))
  file.remove("test_hyper.bed", "test_hypo.bed", "test.seg", "test_z_scores.seg")
})

test_that("compute_z_scores in relaxed conditions and write_output", {
  dmr_table_r <- find_dmrs(tumor_toy_table, control_toy_table,
                                        auc_toy_vector, reference_toy_table)
  dmr_table_r$start <- dmr_table_r$start - 10
  dmr_table_r$end <- dmr_table_r$end + 10

  sample_score_r <- compute_z_scores(tumor_toy_table, control_toy_table,
                                     dmr_table_r, reference_toy_table, 1)
  expect_is(sample_score_r, "list")
  expect_length(sample_score_r, 4)
  write_output("test", dmr_table_r, sample_score_r, 0.8)
  expect_true(file.exists("test_hyper.bed"))
  expect_true(file.exists("test_hypo.bed"))
  expect_true(file.exists("test.seg"))
  expect_true(file.exists("test_z_scores.seg"))
  file.remove("test_hyper.bed", "test_hypo.bed", "test.seg", "test_z_scores.seg")
})

context("long test") ########################################################
test_that("TCGA-ESCA works", {
  data_folder <- "/projects/databases/data/TCGA/harmonized/ESCA/"
  skip_if(file.exists("/projects/packages/ROCkerMeth/skip_long_test"))
  skip_if_not(dir.exists(data_folder))
  x <- round(readRDS(file.path(data_folder, "ESCA_DNAm_TP.rds"))*100)
  y <- round(readRDS(file.path(data_folder, "ESCA_DNAm_NT.rds"))*100)
  auc <- readRDS("/projects/packages/ROCkerMeth/data-raw/pancancer_auc.rds")[,"ESCA"]

  dmr_table <- find_dmrs(x, y, auc, illumina450k_hg19[3:4])
  sample_score <- compute_z_scores(x, y, dmr_table, illumina450k_hg19[3:4], 1)
  write_output(dmr_table, sample_score, path="/projects/packages/ROCkerMeth/esca")
  expect_is(sample_score, "list")
  expect_is(dmr_table, "data.frame")
  expect_is(sample_score, "list")
  expect_true(file.exists("/projects/packages/ROCkerMeth/esca.seg"))
  file.remove("/projects/packages/ROCkerMeth/esca_hyper.bed",
              "/projects/packages/ROCkerMeth/esca_hypo.bed",
              "/projects/packages/ROCkerMeth/esca.seg",
              "/projects/packages/ROCkerMeth/esca_z_scores.seg")
})
