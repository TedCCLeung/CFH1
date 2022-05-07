make_df_BF <- function(

){

  df_DE <- get_DEgenes(edgeR_result, log2FC_threshold = log2FC_threshold, FDR_threshold = FDR_threshold, full = TRUE)
  CPM_genes <- filter_CPM(contrast_list)

  contrasts <- names(CPM_genes)

  df_BF <- do.call("rbind", lapply(contrasts, function(contrast){
    df_1 <- df_DE[(df_DE$contrast == contrast) & (df_DE$ID %in% CPM_genes[[contrast]]), ]
    df_2 <- df_DE[(df_DE$contrast == contrast) & !(df_DE$ID %in% CPM_genes[[contrast]]), ]
    df_2$direction <- rep(0, nrow(df_2))
    return(rbind(df_1, df_2))
  }))

  #utils::write.table(df_BF, file = paste0(dir, "df_BF.tsv"), quote = FALSE, sep= "\t", row.names = FALSE)

  return(df_BF)
}






