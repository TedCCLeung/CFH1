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

  ## Volcano plots ---------------------
  plot_volcano(
    plot_result = edgeR_result[c("SINGLEvCOL", "SINGLEvTRIPLE", "SINGLEvDECOY", "SINGLEvPIF3OX")],
    FDR_threshold = 0.05,
    log2FC_threshold = 1
  ) %>%
    ggplot2::ggsave(file = paste0(dir, "/volcano_single.pdf"), height = 3, width = 6)

  plot_volcano(
    plot_result = plot_result[c("COLvSINGLE", "COLvTRIPLE", "COLvDECOY", "COLvPIF3OX")],
    FDR_threshold = 0.05,
    log2FC_threshold = 1
  ) %>%
    ggplot2::ggsave(file = paste0(dir, "/volcano_COL.pdf"), height = 3, width = 6)


  ## Heatmaps ---------------------
  ## Obtain genes that are differentially expressed in any comparison


}
