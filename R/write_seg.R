#'
write_seg <- function(sample_score, reference_table, outfolder){
  sample_score_molten <- reshape2::melt(sample_score)
  idx <- stringr::str_split((sample_score_molten$Var1), "-", simplify = T)
  genome_mapping <- cbind(reference_table[idx[, 1], c("Chromosome","Genomic_Coordinate")],
        reference_table[idx[, 2],"Genomic_Coordinate"])
  df <- cbind(genome_mapping, sample_score_molten[c("value", "Var2")])
  write.table(df, file = file.path(outfolder, "dmr_z_score.seg"), quote = F, sep = "\t",
              col.names = F, row.names = F)
}
