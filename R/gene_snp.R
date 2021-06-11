#' Match genes with SNPs
#'
#' Internal routine to find features that overlap with SNPs
#'
#' @param snp_tbl tbl of SNPs from \code{query_variants}; see package \code{\link[qtl2]{create_variant_query_func}}
#' @param feature_tbl tbl of feature information from \code{query_genes}; see package \code{\link[qtl2]{create_gene_query_func}}
#' @param feature_snp tbl of feature information from \code{\link{get_feature_snp}}
#'
#' @return tbl of genes covering SNPs
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @importFrom dplyr distinct filter mutate select
#' @importFrom rlang .data
#' 
get_gene_snp <- function(snp_tbl, feature_tbl,
                         feature_snp =
                           get_feature_snp(snp_tbl, feature_tbl, 0)) {
  
  if(is.null(feature_snp))
    return(NULL)
  if(!nrow(feature_snp))
    return(NULL)
  
  out <- dplyr::mutate(
    dplyr::select(
      dplyr::mutate(
        dplyr::distinct(
          dplyr::filter(feature_snp, 
                        .data$type == "gene" & !is.na(.data$Name)), 
          .data$SNP, .data$start, .data$stop, .data$strand, .keep_all=TRUE), 
        gene = .data$Name),
      .data$SNP, .data$pos, .data$lod, .data$gene, .data$start, .data$stop, .data$strand),
    in_gene = (.data$start <= .data$pos & .data$pos <= .data$stop))
  
  class(out) <- c("gene_snp", class(out))
  out
  
}
#' Summary of genes overlapping SNPs
#'
#' @param object tbl of feature information from \code{\link{get_feature_snp}}
#' @param ... additional parameters ignored
#'
#' @return tbl of feature summaries by type
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method summary gene_snp
#' @rdname gene_snp
#' @export
#' @importFrom dplyr group_by n summarize ungroup
summary.gene_snp <- function(object, ...) {
  dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(object, .data$gene, .data$start, .data$stop), 
      snp_count = dplyr::n(),
      min_lod = min(.data$lod),
      max_lod = max(.data$lod)))
}
