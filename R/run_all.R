#' Run all ROCker-meth analysis
#'
#' A simple wrapper to perform a complete analysis with rocker meth. It doesn't
#' return anything but saves two files for further analysis and evaluation.
#'
#' @inheritParams compute_AUC
#' @export
run_all <- function(beta_table, state, reference_table, nclust = 1, NA_thr = 0.5){
  # compute_auc
  auc <- calc_AUC(beta_table, state, nclust, NA_thr)
  # find DMR
  all_dmr <- all_chr_segmentator(auc, beta_table, state, reference_table)
  # compute sample-score
  sample_score <- calc_sample_score(all_dmr, beta_table, state)
  # write bedfile
  write_bed(all_dmr, "/tmp")
  # write segfile
  write_seg(sample_score, reference_table, "/tmp")
}
