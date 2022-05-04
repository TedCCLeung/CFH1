get_DEgenes <- function(
  edgeR_result,
  log2FC_threshold = 1,
  FDR_threshold = 0.05
){

  res_list <- lapply(edgeR_result, function(x){
    return(x[(x$FDR < FDR_threshold) & (abs(x$logFC) >= log2FC_threshold), "genes"])
  })
  names(res_list) <- names(edgeR_result)
  return(res_list)

}
