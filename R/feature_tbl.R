#' Summary of features
#'
#' Show count min and max of features by type
#'
#' @param object tbl of feature information from \code{\link[qtl2]{create_gene_query_func}}
#' @param ... additional parameters passed to methods
#' @param major if \code{TRUE} (default), only summarize genes and exons
#'
#' @return tbl of feature summaries by type
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method summary feature_tbl
#' @rdname feature_tbl
#' @export
#' @importFrom dplyr filter group_by n summarize ungroup
summary.feature_tbl <- function(object, major=TRUE, ...) {
  if(!nrow(object))
    return(NULL)

  if(major)
    object <- dplyr::filter(object,
                            .data$type %in% c("exon","gene","pseudogene","pseudogenic_exon"))
  dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(object, .data$type),
      count = dplyr::n(),
      min_Mbp = min(.data$start),
      max_Mbp = max(.data$stop)))
}
#' Subset of features
#'
#' @param x tbl of feature information from \code{\link{get_feature_snp}}
#' @param start_val,stop_val start and stop positions for subset
#' @param ... additional parameters ignored
#'
#' @return tbl of feature summaries by type
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method subset feature_tbl
#' @rdname feature_tbl
#' @export
#' @importFrom dplyr filter
subset.feature_tbl <- function(x, start_val=0, stop_val=max(x$stop), ...) {
  x <- dplyr::filter(x,
                     .data$start >= start_val,
                     .data$stop <= stop_val)
  class(x) <- unique(c("feature_tbl", class(x)))
  x
}
