#' @importFrom utils write.table
#' 
write_bed <- function(segs, outfolder, fdr_thr = 0.05, zero_based = T){
  out_data <- segs[segs$state != 2 & segs$fdr < 0.05, c("chr", "start", "end", "state", "FDR")]
  if (zero_based){
    out_data[["start"]] <- out_data[["start"]] - 1
    out_data[["end"]] <- out_data[["end"]] - 1
  }
  write.table(out_data, file = file.path(outfolder, "dmr.bed"), quote = F,
              sep = "\t", row.names = F, col.names = F)
}


