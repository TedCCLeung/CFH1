plot_volcano <- function(
  edgeR_result,
  FDR_threshold = 0.05,
  log2FC_threshold = 1
){

  ## Transform the columns and add the CPM information
  edgeR_result_ <- edgeR_result %>%
    lapply(function(x){x$FDR <- -1*log10(x$FDR); x$PValue <- -1*log10(x$PValue); return(x)})

  ## Combine the result
  df_plot <- Reduce(rbind, edgeR_result_)
  df_plot$contrast <- rep(names(edgeR_result_), times = sapply(edgeR_result_, nrow))

  df_plot$DE <- ifelse(
    (df_plot$FDR >= -1*log10(FDR_threshold)) & (abs(df_plot$logFC) >= log2FC_threshold),
    "PASS",
    "FAIL"
  )

  plot <- ggplot2::ggplot(df_plot, ggplot2::aes_string(x = "logFC", y = "FDR")) +
    ggplot2::geom_point(mapping = ggplot2::aes_string(fill = "DE", color = "DE"), size = 0.5, shape = 20, alpha = 0.8) +
    ggplot2::facet_wrap(~contrast, nrow = 1) +
    ggplot2::ylab("-log10(FDR)") +
    #ggplot2::scale_fill_manual(values = c("FAIL" = "#FFFFFF", "PASS" = "#FFFFFF")) +
    ggplot2::scale_color_manual(values = c("FAIL" = "#000000", "PASS" = "#FF0000")) +
    ggplot2::theme(
      legend.position="none",
      #strip.text = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0.0, 0.0, 0.0, 0.0), "cm"),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", size = 0.67, linetype = 1, lineend = "square"),
      axis.text = ggplot2::element_text(size = 7, colour = "black", face = "bold"),
      axis.title = ggplot2::element_text(size = 7, face = "bold")
    )
  return(plot)
}
