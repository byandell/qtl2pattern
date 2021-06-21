#' Merge SNP lod peaks with SNP feature information
#'
#' Merge all SNPs in small region with LOD peaks across multiple phenotype.
#'
#' @param top_snps_tbl tbl from \code{\link{top_snps_all}} or \code{\link[qtl2]{top_snps}}
#' @param snpinfo SNP information table
#' @param out_lmm_snps tbl from \code{\link[qtl2]{scan1}} on SNPs
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
#' example(top_snps_all)
#' 
#' # Download Gene info for DOex from web via RDS
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "c2_genes.rds"), tmpfile, quiet=TRUE)
#' gene_tbl <- readRDS(tmpfile)
#' unlink(tmpfile)
#' 
#' out <- merge_feature(top_snps_tbl, snpinfo, scan_snppr, gene_exon = gene_tbl)
#' summary(out, "pattern")
#'
#' @export
#' @importFrom dplyr arrange distinct filter mutate select
#' @importFrom rlang .data
#' @importFrom qtl2 top_snps
#'
merge_feature <- function(top_snps_tbl, snpinfo, out_lmm_snps, drop=1.5,
                          dropchar=0,
                          gene_exon = get_gene_exon_snp(top_snps_tbl)) {
  phename <- dimnames(out_lmm_snps)[[2]]
  
  if(is.null(gene_exon))
    return(NULL)
  
  ## Add lod by phename to top_snps_tbl
  top_snps_tbl <- dplyr::arrange(
    dplyr::select(
      dplyr::distinct(top_snps_tbl, .data$snp_id, .keep_all=TRUE),
      -.data$pheno),
    .data$pos)
  
  ## Add columns for exons.
  tmp <- top_snps_tbl$pos
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
        qtl2::top_snps(
          subset(out_lmm_snps, lodcolumn = match(i, phename)),
          snpinfo,
          drop=drop),
        .data$snp_id %in% near_snp_id),
      .data$index, .keep_all=TRUE)
    top_snps_tbl[[i]] <- tmp2$lod[match(top_snps_tbl$index, tmp2$index)]
  }
  if(!("consequence" %in% names(top_snps_tbl))) {
    top_snps_tbl$consequence <- "unknown"
  }
  out <- dplyr::select(
    dplyr::mutate(top_snps_tbl,
                  snp_type = abbreviate(.data$consequence,15)),
    -.data$lod)
  
  # haplos
  haplos <- snpinfo_to_haplos(snpinfo)
  attr(out, "haplos") <- haplos
  
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
  haplos <- attr(object, "haplos")
  sum_type <- match.arg(sum_type)
  switch(sum_type,
         "SNP type" = {
           t(c(table(object$snp_type)))
         },
         "pattern" = {
           t(c(table(sdp_to_pattern(object$sdp, haplos))))
         })
}
