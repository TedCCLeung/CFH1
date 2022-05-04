run_heatmap <- function(
  mat,
  genes,
  k = 1,
  distance = "euclidean",
  method = "complete"
){
  mat <- mat[rownames(mat) %in% genes, ]

  ht <- ComplexHeatmap::Heatmap(
    mat,
    clustering_distance_rows = distance,
    clustering_method_rows = method,
    row_km = k,
    show_row_names = FALSE,
    col = circlize::colorRamp2(c(-2, 0, 2), c("#377EB8", "#FFFFFF", "#E41A1C"))
  )

  return(ht)
}
