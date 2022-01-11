context("long test") ########################################################
test_that("TCGA-ESCA works", {
    data_folder <- "/projects/databases/data/TCGA/harmonized/ESCA/"
    skip_if(file.exists("/projects/packages/Rockermeth/skip_long_test"))
    skip_if_not(dir.exists(data_folder))
    x <- round(readRDS(file.path(data_folder, "ESCA_DNAm_TP.rds"))*100)
    y <- round(readRDS(file.path(data_folder, "ESCA_DNAm_NT.rds"))*100)
    auc <- readRDS("/projects/packages/Rockermeth/data-raw/pancancer_auc.rds")[,"ESCA"]

    dmr_table <- suppressMessages(find_dmrs(x, y, auc, illumina450k_hg19, ncores=10))
    sample_score <- suppressMessages(compute_z_scores(x, y, dmr_table, illumina450k_hg19))
    write_dmr(dmr_table, path="/projects/packages/Rockermeth/esca")
    write_z_scores(sample_score, path="/projects/packages/Rockermeth/esca")
    expect_true(file.exists("/projects/packages/Rockermeth/esca.seg"))
    expect_true(file.exists("/projects/packages/Rockermeth/esca_z_scores.seg"))
    file.remove(Sys.glob("/projects/packages/Rockermeth/esca*"))
})
