
DE_discrete_matrix <- function(
  conditions,
  log2FC_threshold = 1,
  FDR_threshold = 0.05
){

  df_DE <- get_DEgenes(edgeR_result, log2FC_threshold = log2FC_threshold, FDR_threshold = FDR_threshold, full = TRUE)
  DE_genes <- get_DEgenes(edgeR_result, log2FC_threshold = log2FC_threshold, FDR_threshold = FDR_threshold)

  df_DE_wide <- tidyr::pivot_wider(df_DE, names_from = "contrast", values_from = "direction")
  df_DE_wide <- df_DE_wide[df_DE_wide$ID %in% unlist(unique(DE_genes)), c("ID", conditions)] %>%
    as.matrix() %>%
    tidyr::replace_na(0)

  return(df_DE_wide)
}
