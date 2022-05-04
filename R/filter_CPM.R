## Function to get genes that are expressed to a certain CPM in both conditions
## Useful because edgeR calculates logFC and p-value even if a gene is only expressed in one condition

filter_CPM <- function(
  contrast_list,
  CPM_threshold = 0
){
  df_CPM <- get_read_matrix(gene_count_matrix, CPM_or_TMM = "CPM")
  CPM_genes <- lapply(contrast_list, function(contrast){
    sample_names <- contrast %>% strsplit(split = "_") %>% sapply(utils::head, 1) %>% unique()
    CPM_filter <- df_CPM[, sample_names] %>% apply(1, function(x){all(x > CPM_threshold)})
    return(df_CPM[CPM_filter, "ID"])
  })
  names(CPM_genes) <- names(contrast_list)
  return(CPM_genes)
}
