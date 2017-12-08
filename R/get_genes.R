#' Get gene in region
#' 
#' Get gene using \code{query_genes}; see \code{\link[qtl2]{create_gene_query_func}}.
#' 
#' @param chr_id chromosome identifier
#' @param start start position in Mbp
#' @param stop  stop position in Mbp
#' @param gene_tbl table of genes from user supplied \code{query_genes}; see \code{\link[qtl2]{create_gene_query_func}}
#' 
#' @export
#' @importFrom dplyr filter
#' 
get_genes <- function(chr_id, start, stop,
                      gene_tbl = query_genes(chr_id, start, stop)) {
  out <- dplyr::filter(
    gene_tbl,
    !is.na(Name))
  class(out) <- c("feature_tbl", class(out))
  out
}