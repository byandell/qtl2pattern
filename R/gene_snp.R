#' Match genes with SNPs
#'
#' Find features that overlap with SNPs
#'
#' @param snp_tbl tbl of SNPs from \code{assoc.map}
#' @param feature_tbl tbl of feature information from \code{query_variants}; see package \code{qtl2db}
#' @param feature_snp tbl of feature information from \code{\link{get_feature_snp}}
#'
#' @return tbl of genes covering SNPs
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{get_gene_snp(snp_tbl, feature_tbl)}
#'
#' @export
#' @importFrom dplyr distinct filter mutate select
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
                        type=="gene" & !is.na(Name)), 
          SNP, start, stop, strand, .keep_all=TRUE), 
        gene=Name),
      SNP, pos, lod, gene, start, stop, strand),
    in_gene = (start<=pos & pos<=stop))
  
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
#' @importFrom dplyr group_by summarize ungroup
summary.gene_snp <- function(object, ...) {
  dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(object, gene, start, stop), 
      snp_count = n(),
      min_lod = min(lod),
      max_lod = max(lod)))
}
