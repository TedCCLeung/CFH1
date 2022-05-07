get_log2FC_mat <- function(

){

  df_TMM <- get_read_matrix(gene_count_matrix, CPM_or_TMM = "TMM")

  mat <- data.frame(
    SINGLEvsCOL = log2(df_TMM$cfh1/df_TMM$COL),
    TRIPLEvsCOL = log2(df_TMM$cfh123/df_TMM$COL),
    DECOYvsCOL = log2(df_TMM$CFHdF/df_TMM$COL),
    PIF3OXvsCOL = log2(df_TMM$PIF3OX/df_TMM$COL)
  ) %>% as.matrix() %>%
    magrittr::set_rownames(df_TMM$ID) %>%
    remove_inf_NA_rows()

  return(mat)
}
