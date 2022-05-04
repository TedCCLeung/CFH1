get_read_matrix <- function(
  gene_count_matrix,
  CPM_or_TMM = "TMM",
  mean = TRUE
){

  ## Select the data
  gene_count_matrix <- gene_count_matrix[order(gene_count_matrix$gene_id), ]
  gene_names <- gene_count_matrix$gene_id %>% strsplit(split = "\\|") %>% sapply(utils::head, 1)
  geneCounts <- as.matrix(gene_count_matrix[, 2:ncol(gene_count_matrix)])
  contrast <- colnames(geneCounts)
  design_matrix <- data.frame(
    Sample = contrast %>% as.factor(),
    Genotype = contrast %>% strsplit(split = "_") %>% lapply(utils::head, 1) %>% unlist() %>% as.factor()
  )
  design <- stats::model.matrix(~design_matrix$Genotype)
  rownames(design) <- design_matrix$Sample

  dge <- edgeR::DGEList(counts = geneCounts, genes = gene_names, group = design_matrix$Genotype)
  #dge <- dge[edgeR::filterByExpr(dge), , keep.lib.sizes = FALSE]

  if (!mean){
    if (CPM_or_TMM == "TMM"){
      dge <- edgeR::calcNormFactors(dge, method = 'TMM')
      return(edgeR::cpm(dge))
    } else {
      return(edgeR::cpm(dge, normalized.lib.sizes = FALSE, prior.count = 0))
    }
  }


  if (mean){
    if (CPM_or_TMM == "TMM"){
      dge <- edgeR::calcNormFactors(dge, method = 'TMM')
      df_count <- edgeR::cpm(dge)
      conds <- unique(design_matrix$Genotype) %>% as.character()
      df_mean <- lapply(conds, function(x){df_count[, design_matrix$Genotype == x] %>% apply(1, mean)}) %>%
        as.data.frame()
      colnames(df_mean) <- conds
      df_mean$ID <- dge$genes[, 1] %>% as.character()
      df_mean <- df_mean[, c("ID", colnames(df_mean)[1:ncol(df_mean)-1])]
      return(df_mean)
    } else {
      df_count <- edgeR::cpm(dge, normalized.lib.sizes = FALSE, prior.count = 0)
      conds <- unique(design_matrix$Genotype) %>% as.character()
      df_mean <- lapply(conds, function(x){df_count[, design_matrix$Genotype == x] %>% apply(1, mean)}) %>%
        as.data.frame()
      colnames(df_mean) <- conds
      df_mean$ID <- dge$genes[, 1] %>% as.character()
      df_mean <- df_mean[, c("ID", colnames(df_mean)[1:ncol(df_mean)-1])]
      return(df_mean)
    }
  }

}
