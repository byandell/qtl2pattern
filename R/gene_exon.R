#' Get exons for set of genes
#'
#' Match up exon start,stop,strand with genes. Use \code{query_genes} to find features; see \code{\link[qtl2]{create_gene_query_func}}.
#'
#' @param top_snps_tbl table from \code{\link[qtl2]{top_snps}}
#' @param feature_tbl table of features from \code{query_genes}; see \code{\link[qtl2]{create_gene_query_func}}
#'
#' @return tbl of exon and gene features
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{get_gene_exon_snp(top_snps_tbl)}
#'
#' @export
#' @rdname gene_exon
#' @importFrom dplyr arrange desc distinct everything mutate select
get_gene_exon_snp <- function(top_snps_tbl,
                              feature_tbl = query_genes(chr_id, range_Mbp[1], range_Mbp[2])) {
  ## Only need distinct snp_id.
  top_snps_tbl <- dplyr::arrange(
    dplyr::select(
      dplyr::distinct(top_snps_tbl, snp_id, .keep_all=TRUE),
      -pheno),
    pos)

  if(is.null(top_snps_tbl))
    return(NULL)
  if(!nrow(top_snps_tbl))
    return(NULL)

  chr_id <- as.character(unique(top_snps_tbl$chr))
  if(length(chr_id) != 1)
    stop("need exactly 1 chromosome in top_snps_tbl")
  range_Mbp <- range(top_snps_tbl$pos) + c(-1,1) * 0.005
  gene_snp <- get_gene_snp(
    dplyr::select(
      top_snps_tbl, 
      snp_id,pos,lod),
    feature_tbl)
  get_gene_exon(feature_tbl, gene_snp)
}

#' Get exons for set of genes
#'
#' Match up exon start,stop,strand with genes.
#'
#' @param feature_tbl tbl of features from \code{query_variants}; see \code{\link[qtl2]{create_variant_query_func}}
#' @param gene_snp tbl of genes with SNPs IDs from \code{\link{match_feature_snp}}
#'
#' @return tbl of exon and gene features
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{get_gene_exon(feature_snp)}
#'
#' @export
#' @rdname gene_exon
#' @importFrom dplyr bind_rows distinct filter
get_gene_exon <- function(feature_tbl, gene_snp) {
  
  if(is.null(gene_snp))
    return(NULL)
  
  ## Need to get unique genes -- duplication with SNPs.
  gene_snp <- dplyr::distinct(gene_snp, gene, .keep_all=TRUE)
  
  if(!nrow(gene_snp)) {
    return(NULL)
  }

  ## Use gene name as Name from feature_tbl.
  exons <- list()
  for(exoni in seq_len(nrow(gene_snp))) {
    genei <- gene_snp$gene[exoni]
    ## get genei and exons spanning genei
    exons[[genei]] <- dplyr::filter(
      dplyr::filter(feature_tbl,
                    type %in% c("exon","gene"),
                    start >= gene_snp$start[exoni],
                    stop <= gene_snp$stop[exoni]),
      (!is.na(Name) & Name==genei) | type=="exon")
    strandi <- gene_snp$strand[exoni]
    if(strandi != "." & !is.na(strandi))
      exons[[genei]] <- dplyr::filter(exons[[genei]], strand==strandi)
  }
  out <- dplyr::distinct(
    dplyr::bind_rows(exons, .id="gene"),
    start, stop, strand, .keep_all=TRUE)
  class(out) <- c("gene_exon", class(out))

  out
}

#' Summary of exons for a gene with SNPs
#'
#' Returns table of gene and its exons.
#'
#' @param gene_exon tbl of feature information from \code{\link{get_gene_exon}}
#' @param gene_name name of gene as character string
#' @param top_snps_tbl table of top SNPs in region from \code{\link[qtl2]{top_snps}}
#' @param extra extra region beyond gene for SNPs (in Mbp)
#'
#' @return tbl of summary
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method summary gene_exon
#' @rdname gene_exon
#' @export
#' @importFrom dplyr arrange desc distinct filter group_by mutate n select summarize ungroup
summary.gene_exon <- function(gene_exon, gene_name=NULL,
                              top_snps_tbl = NULL,
                              extra = 0.005) {
  ## Want to add columns for each phenotype
  ## with number of SNPs within extra of each gene.
  ## How to do this cleverly?
  if(!is.null(gene_name)) {
    out <- dplyr::distinct(
      dplyr::arrange(
        dplyr::select(
          dplyr::filter(gene_exon, gene==gene_name),
          gene, source, type, start, stop, strand),
        dplyr::desc(type)),
      start, stop, strand, .keep_all=TRUE)
  } else {
    out <- dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(
          dplyr::distinct(
            dplyr::filter(gene_exon, type != "gene"),
            start, stop, strand, .keep_all=TRUE),
          gene),
        exons = dplyr::n(),
        min.len = min(stop-start),
        max.len = max(stop-start),
        sum.len = sum(stop-start),
        min_Mbp = min(start),
        max_Mbp = max(stop),
        strand = strand[1]))
    out <- dplyr::select(
      out,
      gene, exons, strand, min_Mbp, max_Mbp, dplyr::everything())
    
    if(!is.null(top_snps_tbl)) {
      ## Goal: add columns to out for each pheno in top_snps_tbl.
      ## Column should have number of SNPs within extra of gene.
      top_snps_tbl <- dplyr::select(top_snps_tbl, pheno, pos, lod)
      pheno_names <- sort(unique(top_snps_tbl$pheno))
      outlim <- out[,c("min_Mbp","max_Mbp")]
      outlim[,1] <- outlim[,1] - extra
      outlim[,2] <- outlim[,2] + extra
      for(pheno_val in pheno_names) {
        out[[pheno_val]] <-
          apply(outlim, 1,
                function(x,y) {
                  in_region <- y$pos >= x[1] & y$pos <= x[2]
                  if(any(in_region))
                    max(y$lod[in_region])
                  else
                    NA
                },
                dplyr::filter(top_snps_tbl, pheno == pheno_val))
      }
      tmp <- out[pheno_names]
      out$max_lod <- apply(tmp, 1, max)
      out <- dplyr::arrange(
        dplyr::select(
          out,
          gene, max_lod, dplyr::everything()),
        dplyr::desc(max_lod))
    }
  }
  out
}
#' @rdname gene_exon
#' @export
subset.gene_exon <- function(x, gene_val, ...) {
  x <- dplyr::filter(x, gene %in% gene_val)
  class(x) <- c("gene_exon", class(x))
  x
}
