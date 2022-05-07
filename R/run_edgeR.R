#' Perform KEGG term enrichment
#'
#' @importFrom magrittr %>%
#'
#' @param contrast Character vector. Column names from the gene_count_matrix.
#' @param prefix Character. Add prefix to some of the column names in the output
#'
#' @return Data frame of edgeR result
#' @export


run_edgeR <- function(
  contrast
){

  ## Select the data
  gene_names <- gene_count_matrix$gene_id %>% strsplit(split = "\\|") %>% sapply(utils::head, 1)
  geneCounts <- as.matrix(gene_count_matrix[, contrast])
  design_matrix <- data.frame(
    Sample = contrast,
    Genotype = contrast %>% strsplit(split = "_") %>% lapply(utils::head, 1) %>% unlist() %>% factor(
      levels = contrast %>% strsplit(split = "_") %>% lapply(utils::head, 1) %>% unlist() %>% unique())
  )

  design <- stats::model.matrix(~design_matrix$Genotype)
  rownames(design) <- design_matrix$Sample

  dge <- edgeR::DGEList(counts = geneCounts, genes = gene_names, group = design_matrix$Genotype)
  dge <- dge[edgeR::filterByExpr(dge), , keep.lib.sizes = FALSE]
  ## Normalize to library size
  dge <- edgeR::calcNormFactors(dge)
  ## Dispersion
  dge <- edgeR::estimateDisp(dge, design, robust=TRUE)

  #edgeR::plotBCV(dge)
  #edgeR::plotMDS.DGEList(dge)

  fit <- edgeR::glmFit(dge, design)
  lrt <- edgeR::glmLRT(fit)
  edgeR_result <- edgeR::topTags(lrt, n = Inf)$table[, c("genes", "PValue", "FDR", "logFC")]
  colnames(edgeR_result) <- c("ID", "PValue", "FDR", "logFC")
  edgeR_result$logFC <- edgeR_result$logFC*-1
  return(edgeR_result)
}



