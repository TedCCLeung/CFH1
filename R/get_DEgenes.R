get_DEgenes <- function(
  edgeR_result,
  log2FC_threshold = 1,
  FDR_threshold = 0.05,
  full = FALSE
){

  if (full){

    up_list <- lapply(edgeR_result, function(x){
      return(x[(x$FDR < FDR_threshold) & (x$logFC >= log2FC_threshold), "ID"])
    })

    neutral_list <- lapply(edgeR_result, function(x){
      return(x[(x$FDR >= FDR_threshold) | (abs(x$logFC) < log2FC_threshold), "ID"])
    })

    down_list <- lapply(edgeR_result, function(x){
      return(x[(x$FDR < FDR_threshold) & (x$logFC <= -1*log2FC_threshold), "ID"])
    })

    df_up <- data.frame(
      ID = unlist(up_list),
      contrast = rep(names(edgeR_result), times = lengths(up_list)),
      direction = rep(1, length(unlist(up_list)))
    )

    df_down <- data.frame(
      ID = unlist(down_list),
      contrast = rep(names(edgeR_result), times = lengths(down_list)),
      direction = rep(-1, length(unlist(down_list)))
    )

    df_neutral <- data.frame(
      ID = unlist(neutral_list),
      contrast = rep(names(edgeR_result), times = lengths(neutral_list)),
      direction = rep(0, length(unlist(neutral_list)))
    )

    df_res <- do.call("rbind", list(df_up, df_neutral, df_down))
    row.names(df_res) <- NULL

    return(df_res)

  } else {

    res_list <- lapply(edgeR_result, function(x){
      return(x[(x$FDR < FDR_threshold) & (abs(x$logFC) >= log2FC_threshold), "ID"])
    })
    names(res_list) <- names(edgeR_result)
    return(res_list)

  }

}
