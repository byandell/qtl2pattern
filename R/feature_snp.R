#' Match features with SNPs
#'
#' Find features that overlap with SNPs
#'
#' @param snp_tbl tbl of SNPs from \code{assoc.map}
#' @param feature_tbl tbl of feature information from \code{\link[qtl2]{create_gene_query_func}}
#' @param extend extend region for SNPs in Mbp (default 0.005)
#'
#' @return tbl of features covering SNPs
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @importFrom dplyr arrange bind_rows distinct filter group_by inner_join 
#' mutate select summarize ungroup
#' @importFrom rlang .data
#'
get_feature_snp <- function(snp_tbl, feature_tbl, extend=0.005) {
  feature_tbl <- dplyr::filter(feature_tbl, 
                               .data$type %in% c("exon","gene","pseudogene","pseudogenic_exon"))
  if(!nrow(feature_tbl))
    return(NULL)

  ## Pull features that overlap with SNPs.
  ## Group features by Dbxref for each gene and get bp range.
  genes <- dplyr::mutate(
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(feature_tbl, .data$Dbxref), 
        min_bp = min(.data$start) - extend,
        max_bp = max(.data$stop) + extend)),
    # Catch any missing values for Dbxref
    Dbxref = ifelse(is.na(.data$Dbxref), "NA", .data$Dbxref))
  if(!nrow(genes))
    return(NULL)
  
  ## For each gene, get rows of snp_tbl in bp range.
  tmpfn <- function(x,snp_tbl) {
    keep <- (x[1] <= snp_tbl$pos) & (x[2] >= snp_tbl$pos)
    snp_tbl[keep,]
  }
  tmp <- apply(as.matrix(genes[,2:3]), 1, tmpfn, snp_tbl)
  names(tmp) <- genes$Dbxref
  tmp <- dplyr::bind_rows(tmp, .id = "Dbxref")
  out <- 
    dplyr::arrange(
      dplyr::select(
        dplyr::mutate(
          dplyr::inner_join(
            dplyr::distinct(
              feature_tbl,
              .data$start, .data$stop, .data$strand, .data$type, .keep_all=TRUE),
            tmp, by = "Dbxref"),
          SNP = .data$snp_id), 
        -.data$snp_id),
      .data$Name, .data$type)

  if(!nrow(out))
    return(NULL)

  class(out) <- c("feature_snp", class(out))
  out
}

#' Summary of features with SNP information
#'
#' @param object tbl of feature information from \code{\link{get_feature_snp}}
#' @param ... additional parameters ignored
#'
#' @return tbl of feature summaries by type
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method summary feature_snp
#' @rdname feature_snp
#' @export
#' @importFrom dplyr group_by n summarize ungroup
#' 
summary.feature_snp <- function(object, ...) {
  dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(object, .data$type), 
      count = dplyr::n(),
      minbp = min(.data$start),
      maxbp = max(.data$stop)))
}

