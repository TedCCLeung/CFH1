remove_inf_NA_rows <- function(x){
  inf_filter <- apply(x, 1, function(y){(!is.infinite(sum(y))) & (!is.na(sum(y)))})
  return(x[inf_filter, ])
}
