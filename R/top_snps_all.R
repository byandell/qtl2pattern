#' Top SNPs for all phenotypes scanned
#'
#' Separate fine mapping scans by allele pattern.
#'
#' @param scan1_output output of linear mixed model for \code{phename} (see \code{\link[qtl2scan]{scan1}})
#' @param drop include all SNPs within \code{drop} of max LOD (default 1.5)
#' @param show_all_snps show all SNPs if \code{TRUE}
#'
#' @return table of top_snps at maximum lod for \code{pattern}
#'
#' @export
#' @importFrom dplyr filter inner_join select
#' @importFrom tidyr gather
#' @importFrom CCSanger convert_bp
#'
top_snps_all <- function (scan1_output, drop = 1.5, show_all_snps = TRUE)
{
    map <- scan1_output$map
    if (is.null(map))
        stop("No map found")
    snpinfo <- scan1_output$snpinfo
    if (is.null(snpinfo))
        stop("No snpinfo found")
    chr <- names(map)
    if (length(chr) > 1) {
        warning("Considering only chromosome ", chr)
        chr <- chr[1]
    }

    ## Following is generalized from qtl2scan::top_snps()
    lod_df <- as.data.frame(subset(scan1_output, chr = chr)$lod)
    lod_df$index <- seq(nrow(lod_df))
    lod_df$snp_id <- rownames(lod_df)
    lod_df <- tidyr::gather(lod_df, pheno, lod, -snp_id, -index)
    lod_df <- dplyr::filter(lod_df, lod > max(lod, na.rm = TRUE) - drop)

    snpinfo <- snpinfo[[chr]]

    if (show_all_snps) {
      snpinfo <- dplyr::inner_join(snpinfo,
                                   dplyr::select(lod_df, -snp_id),
                                   by = "index")
    }
    else {
      snpinfo <- dplyr::inner_join(snpinfo,
                                   dplyr::select(lod_df, -index),
                                   by = "snp_id")
    }
    class(snpinfo) <- c("top_snps_all", class(snpinfo))
    snpinfo
}

#' Summary of top SNP object across phenotypes
#'
#' @param object object of class \code{top_snps_tbl}
#' @param sum_type type of summary (one of "range","peak","best")
#' @param shorten_char number of characters to shorten pheno name
#' @param ... other arguments not used
#'
#' @return table summary
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{summary(object)}
#'
#' @method summary top_snps_all
#' @rdname top_snps_all
#' @export
#' @importFrom dplyr arrange desc filter group_by mutate select summarize ungroup
#' @importFrom CCSanger sdp_to_pattern
#'
summary.top_snps_all <- function(object, sum_type=c("range","peak","best"),
                                 shorten_char=0,
                                 ...) {
  sum_type <- match.arg(sum_type)
  switch(sum_type,
         best = { ## Top SNPs across all phenotypes.
           dplyr::mutate(
             dplyr::arrange(
               dplyr::select(object, -index,-sdp),
               dplyr::desc(lod)),
             pheno = shorten_phename(pheno, shorten_char),
             csq = abbreviate(csq, 15))
         },
         range = {
           ## Most frequent SNP patterns within 1.5 of max LOD.
           dplyr::arrange(
             dplyr::mutate(
               dplyr::ungroup(
                 dplyr::summarize(
                   dplyr::group_by(object, sdp, pheno),
                   pct = round(100 * n() / nrow(object), 2),
                   min_lod = min(lod),
                   max_lod = max(lod),
                   max_pos = pos_Mbp[which.max(lod)][1],
                   max_snp = snp[which.max(lod)][1])),
               pattern = CCSanger::sdp_to_pattern(sdp)),
             dplyr::desc(pct))
         },
         peak = {
           dplyr::select(
             dplyr::arrange(
               dplyr::mutate(
                 dplyr::ungroup(
                   dplyr::summarize(
                     dplyr::group_by(
                       dplyr::filter(
                         dplyr::group_by(object, chr, sdp),
                         lod == max(lod)),
                       chr, sdp, index, lod),
                     min_Mbp = min(pos_Mbp),
                     max_Mbp = max(pos_Mbp),
                     phenos = paste(shorten_phename(unique(pheno),
                                                    shorten_char),
                                    collapse=","))),
                 pattern = CCSanger::sdp_to_pattern(sdp)),
               dplyr::desc(lod)),
             -chr)
         })
}
#' Subset of features
#'
#' @param x tbl of feature information from \code{\link{get_feature_tbl}}
#' @param start_val,stop_val start and stop positions for subset
#' @param pheno phenotype name(s) for subset
#' @param ... additional parameters ignored
#'
#' @return subset of \code{x}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method subset top_snps_all
#' @rdname top_snps_all
#' @export
#' @importFrom dplyr filter
subset.top_snps_all <- function(x, start_val=0, stop_val=max(x$pos_Mbp),
                                pheno = NULL, ...) {
  x <- dplyr::filter(x,
                     pos_Mbp >= CCSanger::convert_bp(start_val, FALSE),
                     pos_Mbp <= CCSanger::convert_bp(stop_val, FALSE))
  pheno_val <- pheno # need to be different from column name in x
  if(!is.null(pheno_val))
    x <-dplyr::filter(x, pheno %in% pheno_val)
  class(x) <- unique(c("top_snps_tbl", class(x)))
  x
}

