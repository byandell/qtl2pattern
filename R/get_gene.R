#' Get gene in region
#' 
#' Get gene using \code{query_genes}; see \code{\link[qtl2]{create_gene_query_func}}.
#' 
#' @param chr_id chromosome identifier
#' @param start start position in Mbp
#' @param stop  stop position in Mbp
#' 
#' @export
#' @importFrom dplyr filter
#' 
get_gene <- function(chr_id, start, stop) {
  out <- dplyr::filter(
    query_genes(chr_id, start, stop),
    !is.na(Name))
  class(out) <- c("feature_tbl", class(out))
  out
}