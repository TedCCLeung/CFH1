make_summary <- function(
  edgeR_result,
  contrast_list
){

  df_edgeR <- Map(
    function(x, y){colnames(x)[2:ncol(x)] <- paste0(y, "-", colnames(x)[2:ncol(x)]); return(x)},
    edgeR_result,
    names(contrast_list)
  ) %>%
    purrr::reduce(dplyr::full_join, by = "ID")

  rowSum_result <- data.frame(
    ID = gene_count_matrix$gene_id %>% strsplit(split = "\\|") %>% sapply(utils::head, 1),
    COL = rowSums(gene_count_matrix[,c("COL_1", "COL_2", "COL_3")]),
    SINGLE = rowSums(gene_count_matrix[,c("cfh1_1", "cfh1_2", "cfh1_3")]),
    TRIPLE = rowSums(gene_count_matrix[,c("cfh123_1", "cfh123_2", "cfh123_3")]),
    DECOY = rowSums(gene_count_matrix[,c("CFHdF_1", "CFHdF_2", "CFHdF_3")]),
    PIF3OX = rowSums(gene_count_matrix[,c("PIF3OX_1", "PIF3OX_2", "PIF3OX_3")]),
    TOTAL = rowSums(gene_count_matrix[,2:ncol(gene_count_matrix)])
  )
  rowSum_result <- rowSum_result[rowSums(rowSum_result[, 2:6]) > 0, ]

  annotation_result <- list(get_TAIR("symbols"), get_TAIR("description"), get_TAIR("pubmed")) %>%
    purrr::reduce(dplyr::full_join, by = "ID")

  full_result <- dplyr::right_join(
    annotation_result,
    dplyr::full_join(rowSum_result, df_edgeR, by = "ID"),
    by = "ID"
  )

  return(full_result)
}
