run_master <- function(
  dir = "./data-raw/"
){

  ## Define the contrasts ---------------------
  contrast_list <- list(
    SINGLEvCOL = c("cfh1_1", "cfh1_2", "cfh1_3", "COL_1", "COL_2", "COL_3"),
    SINGLEvTRIPLE = c("cfh1_1", "cfh1_2", "cfh1_3", "cfh123_1", "cfh123_2", "cfh123_3"),
    SINGLEvDECOY = c("cfh1_1", "cfh1_2", "cfh1_3", "CFHdF_1", "CFHdF_2", "CFHdF_3"),
    SINGLEvPIF3OX = c("cfh1_1", "cfh1_2", "cfh1_3", "PIF3OX_1", "PIF3OX_2", "PIF3OX_3"),

    COLvSINGLE = c("COL_1", "COL_2", "COL_3", "cfh1_1", "cfh1_2", "cfh1_3"),
    COLvTRIPLE = c("COL_1", "COL_2", "COL_3", "cfh123_1", "cfh123_2", "cfh123_3"),
    COLvDECOY = c("COL_1", "COL_2", "COL_3", "CFHdF_1", "CFHdF_2", "CFHdF_3"),
    COLvPIF3OX = c("COL_1", "COL_2", "COL_3", "PIF3OX_1", "PIF3OX_2", "PIF3OX_3")
  )

  ## Perform edgeR ---------------------
  edgeR_result <- Map(run_edgeR, contrast_list)

  ## Output a human readable table ---------------------
  df_summary <- make_summary(edgeR_result = edgeR_result, contrast_list = contrast_list)
  utils::write.table(df_summary, file = paste0(dir, "/result.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

  ## Output TMM and CPM table ---------------------
  df_TMM <- get_read_matrix(gene_count_matrix, CPM_or_TMM = "TMM")
  df_CPM <- get_read_matrix(gene_count_matrix, CPM_or_TMM = "CPM")
  utils::write.table(df_TMM, paste0(dir, "TMM_table.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
  utils::write.table(df_CPM, paste0(dir, "CPM_table.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

  ## Get gene lists ---------------------
  DE_genes <- get_DEgenes(edgeR_result, FDR_threshold = 0.05, log2FC_threshold = 1)
  CPM_genes <- filter_CPM(contrast_list)
  ## Get bona fide genes that are really DE because they are present in both comparisons
  BF_genes <- Map(function(x, y){return(intersect(x, y))}, DE_genes, CPM_genes); names(BF_genes) <- names(DE_genes)
  ## filter edgeR result against these genes
  filtered_edgeR_result <- Map(function(x, y){return(x[x$genes %in% CPM_genes[[y]], ])}, edgeR_result, names(edgeR_result))

  ## Volcano plots ---------------------
  plot_volcano(
    edgeR_result = filtered_edgeR_result[c("SINGLEvCOL", "SINGLEvTRIPLE", "SINGLEvDECOY", "SINGLEvPIF3OX")],
    FDR_threshold = 0.05,
    log2FC_threshold = 1
  ) %>%
    ggplot2::ggsave(file = paste0(dir, "/volcano_single.pdf"), height = 3, width = 6)

  plot_volcano(
    edgeR_result = filtered_edgeR_result[c("COLvSINGLE", "COLvTRIPLE", "COLvDECOY", "COLvPIF3OX")],
    FDR_threshold = 0.05,
    log2FC_threshold = 1
  ) %>%
    ggplot2::ggsave(file = paste0(dir, "/volcano_COL.pdf"), height = 3, width = 6)


  ## Heatmaps ---------------------

  mat_vsCOL <- data.frame(
    cfh1vsCOL = log2(df_TMM$cfh1/df_TMM$COL),
    cfh123vsCOL = log2(df_TMM$cfh123/df_TMM$COL),
    CFHdFvsCOL = log2(df_TMM$CFHdF/df_TMM$COL),
    PIF3OXvsCOL = log2(df_TMM$PIF3OX/df_TMM$COL)
  ) %>% as.matrix() %>%
    magrittr::set_rownames(df_TMM$ID) %>%
    remove_inf_NA_rows()

  grDevices::pdf(paste0(dir, "heatmap.pdf"), height = 6, width = 2.5)
  run_heatmap(
    mat_vsCOL,
    genes = unlist(DE_genes[c("COLvSINGLE", "COLvTRIPLE", "COLvPIF3OX", "COLvDECOY")]),
    method = "complete",
    distance = "euclidean"
    ) %>% ComplexHeatmap::draw()
  grDevices::dev.off()


  ## Without DECOY -----
  mat_vsCOL_nodecoy <- data.frame(
    cfh1vsCOL = log2(df_TMM$cfh1/df_TMM$COL),
    cfh123vsCOL = log2(df_TMM$cfh123/df_TMM$COL),
    #CFHdFvsCOL = log2(df_TMM$CFHdF/df_TMM$COL),
    PIF3OXvsCOL = log2(df_TMM$PIF3OX/df_TMM$COL)
  ) %>% as.matrix() %>%
    magrittr::set_rownames(df_TMM$ID) %>%
    remove_inf_NA_rows()

  grDevices::pdf(paste0(dir, "heatmap_noDECOY.pdf"), height = 6, width = 2.5)
  run_heatmap(
    mat_vsCOL_nodecoy,
    genes = unlist(DE_genes[c("COLvSINGLE", "COLvTRIPLE", "COLvPIF3OX")]),
    method = "complete",
    distance = "euclidean"
  ) %>% ComplexHeatmap::draw()
  grDevices::dev.off()

}
