context("AUC functions") #####################################################

test_that("single_AUC returns correct values", {
  a <- as.integer(c(runif(10, 51, 100), rep(NA, 2)))
  b <- as.integer(c(runif(10,  0,  50), rep(NA, 2)))
  d <- as.integer(c(runif(10, 51, 100), rep(NA, 1000)))
  e <- as.integer(c(runif(10,  0,  50), rep(NA, 1000)))
  auc <- single_AUC(c(a, b), c(rep(TRUE, length(a)), rep(FALSE, length(b))), .2)
  expect_equal(auc, 1)
  auc <- single_AUC(c(a, b), c(rep(TRUE, length(a)), rep(FALSE, length(b))), 0)
  expect_equal(auc, NA)
  auc <- single_AUC(c(d, e), c(rep(TRUE, length(d)), rep(FALSE, length(e))), .9)
  expect_equal(auc, NA)

  auc <- single_AUC(c(tumor_toy_table[1,], control_toy_table[1,]), na_threshold=0,
    state = c(rep(TRUE, ncol(tumor_toy_table)), rep(FALSE, ncol(control_toy_table))))
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

test_that("fix_short_segments works", {
  x <- c(3,4)
  expect_error(fix_short_segments(x, 5), "Elements")

  x <- c(3,3,3,3,3,3,2,2,2,1,1,2,3,3,3,2,2,1,1,1,1,1,1)
  y <- c(rep(3,6),rep(2,11),rep(1,6))
  expect_equal(fix_short_segments(x, 5), y)
})

test_that("meth_state_finder works", {
  idx_chr <- which(reference_toy_table[[1]] == "2")
  auc_sd <- sd(auc_toy_vector, na.rm = TRUE)
  x <- auc_toy_vector[idx_chr]
  y <- reference_toy_table[[2]][idx_chr]
  idx_not_NA <- which(!is.na(x))
  meth_states <- meth_state_finder(x[idx_not_NA], y[idx_not_NA], auc_sd, 5,
    pt_start = 0.05, normdist = 1e5, ratiosd = 0.4)
  set.seed(1)
  expect_equal(sample(meth_states, 10), c(2, 2, 2, 2, 1, 2, 2, 1, 2, 2))

})

context("segmentator functions") ######################################
test_that("segmentator returns correct table", {
  idx_chr <- which(reference_toy_table[[1]] == "2")
  auc_sd <- sd(auc_toy_vector, na.rm = TRUE)
  x <- auc_toy_vector[idx_chr]
  y <- reference_toy_table[[2]][idx_chr]
  idx_not_NA <- which(!is.na(x))
  meth_states <- meth_state_finder(x[idx_not_NA], y[idx_not_NA], auc_sd, 5,
    pt_start = 0.05, normdist = 1e5, ratiosd = 0.4)

  tumor_toy_beta_mean <- apply(tumor_toy_table[idx_chr,], 1, mean, na.rm = TRUE)
  normal_beta_mean <- apply(control_toy_table[idx_chr,], 1, mean, na.rm = TRUE)
  single_chr_segs <- segmentator(meth_states, tumor_toy_beta_mean, normal_beta_mean)
  expect_is(single_chr_segs, "data.frame")
  expect_length(single_chr_segs, 4)
  expect_identical(names(single_chr_segs), c("nseg", "state", "avg_beta_diff", "p_value"))
})

test_that("whole_genome_segmentator works", {
  dmr_table <- whole_genome_segmentator(tumor_toy_table, control_toy_table,
    auc_toy_vector, reference_toy_table)
  expect_is(dmr_table, "data.frame")
  expect_length(dmr_table, 8)
  expect_identical(names(dmr_table),
    c("chr", "start", "end", "nseg", "state", "avg_beta_diff", "p_value", "q_value"))
})

context("ouput") ##############################################################
test_that("compute_z_score and write_output works", {
  dmr_table <- whole_genome_segmentator(tumor_toy_table, control_toy_table,
    auc_toy_vector, reference_toy_table)
  sample_score <- compute_z_score(tumor_toy_table, control_toy_table,
    dmr_table, reference_toy_table)
  expect_is(sample_score, "list")
  expect_length(sample_score, 3)
  write_output("test", dmr_table, sample_score, 0.8)
  expect_true(file.exists("test_hyper.bed"))
  expect_true(file.exists("test_hypo.bed"))
  expect_true(file.exists("test.seg"))
  expect_true(file.exists("test_z_scores.seg"))
  file.remove("test_hyper.bed", "test_hypo.bed", "test.seg", "test_z_scores.seg")
  write_output("test", dmr_table, sample_score, 0.8)
})

context("long test") ########################################################
test_that("TCGA-ESCA works", {
  data_folder <- "/projects/databases/data/TCGA/clean/harmonized/ESCA/"
  skip_if_not(dir.exists(data_folder))
  skip("just skip")
  x <- round(readRDS(file.path(data_folder, "ESCA_DNAm_TP.rds"))*100)
  y <- round(readRDS(file.path(data_folder, "ESCA_DNAm_NT.rds"))*100)
  auc <- readRDS("/projects/packages/ROCkerMeth/data-raw/pancancer_auc.rds")[,"ESCA"]
  load("/projects/packages/PAMES/data/illumina450k_hg19.rda")
  idx <- order(illumina450k_hg19$Chromosome, illumina450k_hg19$Genomic_Coordinate)
  x <- x[idx,]
  y <- y[idx,]
  auc <- auc[idx]
  illumina450k_hg19 <- illumina450k_hg19[idx,]

  dmr_table <- whole_genome_segmentator(x, y, auc, illumina450k_hg19[3:4])
  sample_score <- compute_z_score(x, y, dmr_table, illumina450k_hg19[3:4])
  write_output(dmr_table, sample_score, "/projects/packages/ROCkerMeth/esca")
  expect_is(sample_score, "list")
  expect_is(dmr_table, "data.frame")
  expect_is(sample_score, "list")
  expect_true(file.exists("/projects/packages/ROCkerMeth/esca.seg"))
  file.remove("/projects/packages/ROCkerMeth/esca.seg")
})
