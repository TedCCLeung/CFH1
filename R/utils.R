remove_inf_NA_rows <- function(x){
  inf_filter <- apply(x, 1, function(y){(!is.infinite(sum(y))) & (!is.na(sum(y)))})
  return(x[inf_filter, ])
}

addLeadingZeros <- function(
  num_vec
){
  digit_number <- floor(max(log10(num_vec)+1))
  return(stringr::str_pad(as.character(num_vec), digit_number, pad = "0"))
}
