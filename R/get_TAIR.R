#' Function to generate a TxDb object from the GenomicFeatures object for getting gene promoters.
#'
#' @param annotation Either "description", "pubmed", "symbols
#'
#' @return Returns a TxDb object (from the GenomicFeatures package)
#' @export

get_TAIR <- function(
  annotation = "description"
){

  if (annotation == "description"){

    x <- org.At.tair.db::org.At.tairGENENAME
    # Get the TAIR identifiers that are mapped to a gene name
    mapped_tairs <- AnnotationDbi::mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_tairs])
    df <- data.frame(
      ID = mapped_tairs,
      description = lapply(xx, function(x){paste(x, collapse = "; ")}) %>% unlist()
    )
    return(df)

  } else if (annotation == "pubmed"){

    x <- org.At.tair.db::org.At.tairPMID
    # Get the TAIR identifiers that are mapped to a gene name
    mapped_tairs <- AnnotationDbi::mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_tairs])
    df <- data.frame(
      ID = mapped_tairs,
      PMID = lapply(xx, function(x){paste(x, collapse = "; ")}) %>% unlist()
    )
    return(df)

  } else if (annotation == "symbols"){

    x <- org.At.tair.db::org.At.tairSYMBOL
    # Get the TAIR identifiers that are mapped to a gene name
    mapped_tairs <- AnnotationDbi::mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_tairs])
    df <- data.frame(
      ID = mapped_tairs,
      symbol = lapply(xx, function(x){paste(x, collapse = "; ")}) %>% unlist()
    )
    return(df)
    }

}
