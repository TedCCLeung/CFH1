run_heatmap <- function(
  conditions = c("SINGLEvsCOL", "TRIPLEvsCOL", "PIF3OXvsCOL", "DECOYvsCOL"),
  plot_conditions = NULL,
  distance = "pearson",
  method = "complete",
  k,
  return_clusters = FALSE,
  show_row_names = FALSE
){

  if (is.null(plot_conditions)){plot_conditions <- conditions}

  mat_FC <- get_log2FC_mat()[, plot_conditions]
  genes <- DE_genes[conditions] %>% unlist() %>% unique() %>% sort()
  filtered_mat_FC <- mat_FC[rownames(mat_FC) %in% genes, ]

  hclust_out <- filtered_mat_FC %>%
    amap::Dist(method = distance) %>%
    stats::hclust(method = method)

  ht <- ComplexHeatmap::Heatmap(
    filtered_mat_FC,
    cluster_rows = dendextend::color_branches(hclust_out, k = k),
    show_row_names = show_row_names,
    row_split = k,
    col = circlize::colorRamp2(c(-2, 0, 2), c("#377EB8", "#FFFFFF", "#E41A1C")),
    use_raster = TRUE
  )

  if (return_clusters){

    clusterlist = ComplexHeatmap::row_order(ht) %>%
      suppressWarnings() %>%
      suppressMessages()

    df_out <- data.frame(
      ID = rownames(filtered_mat_FC)[unlist(clusterlist)],
      cluster = rep(addLeadingZeros(1:length(clusterlist)), times = lengths(clusterlist))
    )
    return(df_out)

  } else {
    return(ht)
  }

}
