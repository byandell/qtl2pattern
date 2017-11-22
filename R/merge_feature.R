#' Merge SNP lod peaks with SNP feature information
#'
#' Merge all SNPs in small region with LOD peaks across multiple phenotype.
#'
#' @param top_snps_tbl tbl from \code{\link{get_top_snps_tbl}} or \code{\link[qtl2scan]{top_snps}}
#' @param snpinfo SNP information table
#' @param out_lmm_snps tbl from \code{\link[qtl2scan]{scan1}} on SNPs
#' @param drop include LOD scores within \code{drop} of max for each phenotype
#' @param dropchar number of characters to drop on phenames
#' @param gene_exon tbl from \code{\link{get_gene_exon_snp}}
#'
#' @return tbl with added information on genes and exons
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @examples
#' \dontrun{merge_feature(...)}
#'
#' @export
#' @importFrom dplyr arrange distinct filter mutate select
#' @importFrom qtl2scan top_snps
#' @importFrom CCSanger convert_bp get_gene_exon_snp
#'
merge_feature <- function(top_snps_tbl, snpinfo, out_lmm_snps, drop=1.5,
                          dropchar=0,
                          gene_exon =
                            CCSanger::get_gene_exon_snp(top_snps_tbl)) {
  phename <- dimnames(out_lmm_snps)[[2]]

  ## Add lod by phename to top_snps_tbl
  top_snps_tbl <- dplyr::arrange(
    dplyr::select(
      dplyr::distinct(top_snps_tbl, snp_id, .keep_all=TRUE),
      -pheno),
    pos_Mbp)

  ## Add columns for exons.
  tmp <- CCSanger::convert_bp(top_snps_tbl$pos_Mbp)
  ins <- outer(gene_exon$start, tmp, "<=") &
    outer(gene_exon$stop, tmp, ">=")
  ## SNP position should be in 1 (or more if splice variant) exon(s).
  top_snps_tbl$exon_ct <- apply(ins,2,sum)
  top_snps_tbl$exon_id <- apply(ins, 2, function(x) paste(which(x), collapse=";"))

  ## SNP IDs for top_snps_tbl
  near_snp_id <- top_snps_tbl$snp_id

  ## Merge columns for LODs by phenotype.
  for(i in phename) {
    tmp2 <- dplyr::distinct(
      dplyr::filter(
        qtl2scan::top_snps(
          subset(out_lmm_snps, lodcolumn=match(i, phename)),
          snpinfo,
          drop=drop),
        snp_id %in% near_snp_id),
      index, .keep_all=TRUE)
    top_snps_tbl[[i]] <- tmp2$lod[match(top_snps_tbl$index, tmp2$index)]
  }
  out <- dplyr::select(
    dplyr::mutate(top_snps_tbl,
                  snp_type = abbreviate(csq,15)),
    -lod)
  class(out) <- c("merge_feature", class(out))
  out
}
#' Summary of merge_feature object
#'
#' @param object of class \code{merge_feature}
#' @param sum_type one of \code{c("SNP type","pattern")}
#' @param ... other arguments not used
#'
#' @return table summary
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method summary merge_feature
#' @rdname merge_feature
#' @export
summary.merge_feature <- function(object,
                                  sum_type = c("SNP type","pattern"), ...) {
  sum_type <- match.arg(sum_type)
  switch(sum_type,
         "SNP type" = {
           t(table(object$snp_type))
         },
         "pattern" = {
           t(table(sdp_to_pattern(object$sdp)))
         })
}
