#' Plot information about DMRs
#'
#' A way to visualize Differential Methylated Regions showing beta values, 
#' AUC, and Methylation States
#'
plot_DMR <- function(tumor_table, control_table, auc_vector, genome_version = "hg19"){
  tumor_table <- as.matrix(tumor_table)
  control_table <- as.matrix(control_table)
  beta_table <- as.matrix(cbind(tumor_table, control_table))
  if (nrow(tumor_table) != nrow(control_table)){
    stop("tumor_table and control_table must have the same number of rows")
  }
  tumor_is_fraction <- all(tumor_table >= 0, na.rm = T) && all(tumor_table <= 1, na.rm = T)
  control_is_fraction <- all(control_table >= 0, na.rm = T) && all(control_table <= 1, na.rm = T)
  if (tumor_is_fraction && control_is_fraction) {
    tumor_table <- tumor_table*100
    control_table <- control_table*100
    storage.mode(control_table) <- "integer"
    storage.mode(tumor_table) <- "integer"
  } else {
    stop("tumor_table and control_table must contain decimal values")
  }
  sample_classes <- c(rep("Tumor", ncol(tumor_table)), rep("Control", ncol(control_table)))

  gtrack <- GenomeAxisTrack()  # genomic axis

  # user input
  genome_version <- "hg19"
  transcript <- "NM_030903"
  transcript <- NULL
  # mygene <- "OR2W1"
  chr <- "16"
  coord_start <- 28850000
  coord_end <-   28900000

  gene_model <- makeTxDbFromUCSC(genome=genome_version, tablename="refGene", transcript_ids = transcript)

  reference_table = illumina27k_hg19
  idx_locus <- with(reference_table,
    which(Chromosome == chr &
      Genomic_Coordinate >= coord_start &
        Genomic_Coordinate <= coord_end))
  idx_srt <- with(reference_table[idx_locus,], order(Chromosome, Genomic_Coordinate))

  # chromosome representation
  itrack <- IdeogramTrack(genome = genome_version, chromosome = chr)

  # gene region
  grtrack <- GeneRegionTrack(gene_model, genome = genome_version, chromosome = chr,
    fill = "black", min.height = 10, name = "RefGenes")
  grtrack@dp@pars$fill <- "black"
  grtrack@dp@pars$min.height <- 10

  # main track
  main_track <- with(reference_table[idx_locus[idx_srt],], GRanges(sprintf("chr%s", Chromosome),
    IRanges(Genomic_Coordinate, Genomic_Coordinate), "*"))
  genome(main_track) <- genome_version

  # beta values: tumor and normal
  beta_dtrack <- DataTrack(range = main_track,
    data = t(beta_table[idx_locus[idx_srt], ]),
    name = "Beta-Value",
    type = c("a", "confint"),
    showSampleNames = F,
    groups = sample_classes,
    ylim = c(0, 100))
  # auc scores
  auc_gene <- auc_vector[idx_locus[idx_srt]]
  auc_dtrack <- DataTrack(range = main_track,
    data = auc_gene,
    ylim = c(0, 1),
    type = "p",
    showSampleNames = T,
    name = "AUC")
  # DMR values
  dmr_gene <- dmr_table$state[idx_locus[idx_srt]] - 2  # center around 0
  dmr_gene <- c(1,1,3,3)
  dmr_dtrack <- DataTrack(range = main_track,
    name = "DMR",
    data = dmr_gene,
    ylim = c(-1, 1),
    type = "s",
    showSampleNames = F)

  # pdf(file = "/projects/2017_DMR/outs/DMR_example_ORF.draft.pdf", width = 8, height = 6)
  plotTracks(list(itrack, gtrack, grtrack, beta_dtrack, auc_dtrack, dmr_dtrack),
             from = coord_start, to = coord_end)
  # dev.off()
}
