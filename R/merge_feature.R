#' Merge SNP lod peaks with SNP feature information
#'
#' Merge all SNPs in small region with LOD peaks across multiple phenotype.
#'
#' @param top_snps_tbl tbl from \code{\link{get_top_snps_tbl}} or \code{\link[qtl2scan]{top_snps}}
#' @param out_lmm_snps tbl from \code{\link[qtl2scan]{scan1}} on SNPs
#' @param drop include LOD scores within \code{drop} of max for each phenotype
#' @param dropchar number of characters to drop on phenames
#' @param gene_exon tbl from \code{\link{get_gene_exon_snp}}
#' @param sql_filename path to \code{\link{get_mgi_features}}
#' @param datapath path to use with \code{\link{get_mgi_features}}
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
merge_feature <- function(top_snps_tbl, out_lmm_snps, drop=1.5,
                          dropchar=0,
                          gene_exon =
                            CCSanger::get_gene_exon_snp(top_snps_tbl,
                                                        sql_filename),
                          sql_filename = file.path(datapath,
                                                   "mgi_db.sqlite"),
                          datapath) {
  phename <- dimnames(out_lmm_snps$lod)[[2]]
  chr_id <- names(out_lmm_snps$map)
  if(length(chr_id) != 1)
    stop("need exactly 1 chromosome in top_snps_tbl")

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
    tmp2 <- qtl2scan::top_snps(
      dplyr::distinct(
        dplyr::filter(
          subset(out_lmm_snps, lodcolumn=match(i, phename)), drop=drop),
        snp_id %in% near_snp_id),
      index, .keep_all=TRUE)
    top_snps_tbl[shorten_phename(i, dropchar)] <-
      tmp2$lod[match(top_snps_tbl$index, tmp2$index)]
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
#' @param type one of \code{c("SNP type","pattern")}
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
                                  type = c("SNP type","pattern"), ...) {
  type <- match.arg(type)
  switch(type,
         "SNP type" = {
           t(table(object$snp_type))
         },
         "pattern" = {
           t(table(sdp_to_pattern(object$sdp)))
         })
}
#' Plot of merge_feature object
#'
#' @param x of class \code{merge_feature}
#' @param pheno name of phenotype to be plotted
#' @param plot_by element to plot by (one of \code{c("pattern","consequence")})
#' @param ... other arguments not used
#'
#' @return ggplot2 object
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method plot merge_feature
#' @rdname merge_feature
#'
#' @export
#' @importFrom dplyr filter mutate
#' @importFrom ggplot2 aes facet_wrap geom_jitter ggplot ggtitle xlab ylab
#'
plot.merge_feature <- function(x, pheno, plot_by=c("pattern","consequence"), ...) {
  plot_by <- match.arg(plot_by)
  x <- dplyr::filter(
    dplyr::mutate(x, pattern = sdp_to_pattern(sdp), lod = pheno),
    !is.na(lod))
  switch(plot_by,
         pattern = {
           ggplot2::ggplot(x,
                           ggplot2::aes(x=pos_Mbp,y=lod,col=pattern)) +
             ggplot2::geom_jitter() +
             ggplot2::facet_wrap(~snp_type) +
             ggplot2::xlab("Position in Mbp") +
             ggplot2::ylab("LOD") +
             ggplot2::ggtitle(paste("Top SNPs by Consequence for", pheno))
         },
         consequence = {
           ggplot2::ggplot(x,
                           ggplot2::aes(x=pos_Mbp,y=lod,col=snp_type)) +
             ggplot2::geom_jitter() +
             ggplot2::facet_wrap(~pattern) +
             ggplot2::xlab("Position in Mbp") +
             ggplot2::ylab("LOD") +
             ggplot2::ggtitle(paste("Top SNPs by Allele Pattern for", pheno))
         })
}

