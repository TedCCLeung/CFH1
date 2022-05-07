run_master <- function(
  dir = "./data-raw/",
  FDR_threshold = 0.05,
  log2FC_threshold = 0.5
){

  ## Define the contrasts ---------------------
  contrast_list <- list(
    SINGLEvsCOL = c("cfh1_1", "cfh1_2", "cfh1_3", "COL_1", "COL_2", "COL_3"),
    TRIPLEvsCOL = c("cfh123_1", "cfh123_2", "cfh123_3", "COL_1", "COL_2", "COL_3"),
    DECOYvsCOL = c("CFHdF_1", "CFHdF_2", "CFHdF_3", "COL_1", "COL_2", "COL_3"),
    PIF3OXvsCOL = c("PIF3OX_1", "PIF3OX_2", "PIF3OX_3", "COL_1", "COL_2", "COL_3"),

    PIF3OXvsSINGLE = c("PIF3OX_1", "PIF3OX_2", "PIF3OX_3", "cfh1_1", "cfh1_2", "cfh1_3")
  )

  ## Perform edgeR ---------------------
  edgeR_result <- Map(run_edgeR, contrast_list)

  ## Output TMM and CPM table ---------------------
  df_TMM <- get_read_matrix(gene_count_matrix, CPM_or_TMM = "TMM")
  df_CPM <- get_read_matrix(gene_count_matrix, CPM_or_TMM = "CPM")
  utils::write.table(df_TMM, paste0(dir, "TMM_table.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
  utils::write.table(df_CPM, paste0(dir, "CPM_table.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

  ## Get gene lists ---------------------
  DE_genes <- get_DEgenes(edgeR_result, FDR_threshold = FDR_threshold, log2FC_threshold = log2FC_threshold)
  CPM_genes <- filter_CPM(contrast_list)
  ## Get bona fide genes that are really DE because they are present in both comparisons
  BF_genes <- Map(function(x, y){return(intersect(x, y))}, DE_genes, CPM_genes); names(BF_genes) <- names(DE_genes)
  ## filter edgeR result against these genes
  filtered_edgeR_result <- Map(function(x, y){return(x[x$ID %in% CPM_genes[[y]], ])}, edgeR_result, names(edgeR_result))

  ## Output a human readable table ---------------------
  df_summary <- make_summary(edgeR_result = edgeR_result, contrast_list = contrast_list)
  df_DE <- DE_discrete_matrix(
    conditions = names(contrast_list),
    FDR_threshold = FDR_threshold,
    log2FC_threshold = log2FC_threshold
  ) %>% as.data.frame()
  df_full <- dplyr::full_join(df_summary, df_DE, by = "ID")
  utils::write.table(df_full,
                     file = paste0(dir, "/result.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

  ## Volcano plots ---------------------
  plot_volcano(
    edgeR_result = filtered_edgeR_result[c("PIF3OXvsSINGLE")],
    FDR_threshold = FDR_threshold,
    log2FC_threshold = log2FC_threshold
  ) %>%
    ggplot2::ggsave(file = paste0(dir, "/volcano_pif_single.pdf"), height = 3, width = 2)

  ## Heatmaps ---------------------
  ## Heatmap 1
  grDevices::pdf(paste0(dir, "heatmap_noTriple.pdf"), height = 6, width = 2.5)
  run_heatmap(conditions = c("SINGLEvsCOL", "PIF3OXvsCOL", "DECOYvsCOL"), k = 12) %>%
    ComplexHeatmap::draw()
  grDevices::dev.off()
  run_heatmap(conditions = c("SINGLEvsCOL", "PIF3OXvsCOL", "DECOYvsCOL"), k = 12, return_clusters = TRUE) %>%
    utils::write.csv(paste0(dir, "clusters_noTriple.csv"), row.names = FALSE)

  ## Heatmap 2
  grDevices::pdf(paste0(dir, "heatmap_single_PIF.pdf"), height = 6, width = 2.5)
  run_heatmap(conditions = c("SINGLEvsCOL", "PIF3OXvsCOL"), k = 6) %>%
    ComplexHeatmap::draw()
  grDevices::dev.off()
  run_heatmap(conditions = c("SINGLEvsCOL", "PIF3OXvsCOL"), k = 6, return_clusters = TRUE) %>%
    utils::write.csv(paste0(dir, "clusters_single_PIF.csv"), row.names = FALSE)

  ## Heatmap 3
  grDevices::pdf(paste0(dir, "heatmap_noDecoy.pdf"), height = 6, width = 2.5)
  run_heatmap(conditions = c("SINGLEvsCOL", "TRIPLEvsCOL", "PIF3OXvsCOL"), k = 5) %>%
    ComplexHeatmap::draw()
  grDevices::dev.off()
  run_heatmap(conditions = c("SINGLEvsCOL", "TRIPLEvsCOL", "PIF3OXvsCOL"), k = 5, return_clusters = TRUE) %>%
    utils::write.csv(paste0(dir, "clusters_noDecoy.csv"), row.names = FALSE)

  ## Heatmap 4
  grDevices::pdf(paste0(dir, "heatmap_ALL.pdf"), height = 6, width = 2.5)
  run_heatmap(conditions = c("SINGLEvsCOL", "TRIPLEvsCOL", "PIF3OXvsCOL", "DECOYvsCOL"), k = 10) %>%
    ComplexHeatmap::draw()
  grDevices::dev.off()
  run_heatmap(conditions = c("SINGLEvsCOL", "TRIPLEvsCOL", "PIF3OXvsCOL", "DECOYvsCOL"), k = 10, return_clusters = TRUE) %>%
    utils::write.csv(paste0(dir, "clusters_ALL.csv"), row.names = FALSE)

  ## Heatmap 5
  grDevices::pdf(paste0(dir, "heatmap_single_only.pdf"), height = 6, width = 3.5)
  run_heatmap(conditions = c("SINGLEvsCOL"),
              plot_condition = c("SINGLEvsCOL", "PIF3OXvsCOL", "DECOYvsCOL"),
              k = 2,
              show_row_names = TRUE) %>%
    ComplexHeatmap::draw()
  grDevices::dev.off()
}
