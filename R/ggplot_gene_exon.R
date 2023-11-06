#' Plot of exons for a gene with SNPs
#'
#' Uses \code{\link{gene_exon}} to plot genes, exons, mRNA with SNPs.
#'
#' @param object table of feature information from \code{query_genes}; see \code{\link[qtl2]{create_gene_query_func}}
#' @param top_snps_tbl table from \code{\link[qtl2]{top_snps}}
#' @param plot_now plot now if TRUE
#' @param genes Names of genes in \code{object}
#' @param ... arguments passed along to \code{\link{gene_exon}}
#'
#' @return list of ggplots (see \code{\link{gene_exon}})
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @export
#' 
#' @importFrom dplyr group_by summarize ungroup
#' @importFrom ggplot2 ggtitle
#' @importFrom rlang .data
#' 
#' @rdname gene_exon
#' 
ggplot_gene_exon <- function(object, top_snps_tbl=NULL, plot_now=TRUE,
                           genes = unique(object$gene), ...) {
  p <- list()
  ## Reduce to max lod per SNP
  if(!is.null(top_snps_tbl)) {
    top_snps_tbl <- dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(top_snps_tbl, .data$snp_id, .data$pos),
        lod = max(.data$lod)))
  }
  for(genei in genes) {
    p[[genei]] <-
      ggplot_feature_tbl(
        dplyr::mutate(
          dplyr::filter(object, .data$gene == genei),
          gene = NA),
        top_snps_tbl = top_snps_tbl,
        str_rect="",
        ...) +
      ggplot2::ggtitle(genei)
  }
  #*** NONSTANDARD ***
  if(plot_now & length(p)) {
    for(genei in genes)
      print(p[[genei]])
  }
  invisible(p)
}
#' @param object Object of class \code{gene_exon}.
#' @method autoplot gene_exon
#' @export
#' @export autoplot.gene_exon
#' @rdname gene_exon
#' 
#' @importFrom ggplot2 autoplot
#' 
autoplot.gene_exon <- function(object, ...)
  ggplot_gene_exon(object, ...)
