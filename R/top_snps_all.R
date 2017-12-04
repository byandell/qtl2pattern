#' Top SNPs for all phenotypes scanned
#'
#' Separate fine mapping scans by allele pattern.
#'
#' @param scan1_output output of linear mixed model for \code{phename} (see \code{\link[qtl2]{scan1}})
#' 
#' @param snpinfo Data frame with SNP information with the following
#'     columns (the last three are generally derived from with
#'     \code{\link[qtl2]{index_snps}}):
#' \itemize{
#' \item \code{chr} - Character string or factor with chromosome
#' \item \code{pos} - Position (in same units as in the \code{"map"}
#'     attribute in \code{genoprobs}.
#' \item \code{sdp} - Strain distribution pattern: an integer, between
#'     1 and \eqn{2^n - 2} where \eqn{n} is the number of strains, whose
#'     binary encoding indicates the founder genotypes
#' \item \code{snp_id} - Character string with SNP identifier (if
#'     missing, the rownames are used).
#' \item \code{index} - Indices that indicate equivalent
#'     groups of SNPs.
#' \item \code{intervals} - Indexes that indicate which marker
#'     intervals the SNPs reside.
#' \item \code{on_map} - Indicate whether SNP coincides with a marker
#'     in the \code{genoprobs}
#' }
#'
#' @param drop include all SNPs within \code{drop} of max LOD (default 1.5)
#' @param show_all_snps show all SNPs if \code{TRUE}
#'
#' @return table of top_snps at maximum lod for \code{pattern}
#'
#' @export
#' @importFrom dplyr filter group_by inner_join select ungroup
#' @importFrom tidyr gather
#'
top_snps_all <- function (scan1_output, snpinfo, drop = 1.5, show_all_snps = TRUE,
                          haplos)
{
    if (missing(snpinfo) || is.null(snpinfo))
        stop("No snpinfo found")
    map <- snpinfo_to_map(snpinfo)
    
    if(missing(haplos) || is.null(haplos))
      haplos <- snpinfo_to_haplos(snpinfo)

    chr <- names(map)
    if (length(chr) > 1) {
        warning("Considering only chromosome ", chr)
        chr <- chr[1]
    }

    ## Following is generalized from qtl2::top_snps()
    lod_df <- as.data.frame(subset(scan1_output, map, chr = chr))
    lod_df$index <- unique(snpinfo$index)
    lod_df$snp_id <- rownames(lod_df)
    lod_df <- tidyr::gather(lod_df, pheno, lod, -snp_id, -index)
    maxlod <- max(lod_df$lod, na.rm = TRUE) - drop
    lod_df <- dplyr::ungroup(
      dplyr::filter(
        dplyr::group_by(lod_df, pheno),
        lod >= min(max(lod, na.rm = TRUE), maxlod)
      )
    )
#    lod_df <- dplyr::filter(lod_df, lod > max(lod, na.rm = TRUE) - drop)

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
    attr(snpinfo, "haplos") <- haplos
    class(snpinfo) <- c("top_snps_all", class(snpinfo))
    snpinfo
}

#' Summary of top SNP object across phenotypes
#'
#' @param object object of class \code{top_snps_tbl}
#' @param sum_type type of summary (one of "range","peak","best")
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
#'
summary.top_snps_all <- function(object, sum_type=c("range","peak","best"),
                                 ...) {
  haplos <- attr(object, "haplos")
  sum_type <- match.arg(sum_type)
  switch(sum_type,
         best = { ## Top SNPs across all phenotypes.
           out <- 
             dplyr::arrange(
               dplyr::select(object, -index,-sdp),
               dplyr::desc(lod))
           if(!is.null(out$consequence)) {
             out <- 
               dplyr::mutate(out,
                 consequence = abbreviate(consequence, 15))
           }
           out
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
                   max_pos = pos[which.max(lod)][1],
                   max_snp = snp_id[which.max(lod)][1])),
               pattern = sdp_to_pattern(sdp, haplos)),
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
                     min_Mbp = min(pos),
                     max_Mbp = max(pos),
                     phenos = paste(unique(pheno), collapse=","))),
                 pattern = sdp_to_pattern(sdp, haplos)),
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
subset.top_snps_all <- function(x, start_val=0, stop_val=max(x$pos),
                                pheno = NULL, ...) {
  haplos <- attr(x, "haplos")
  x <- dplyr::filter(x,
                     pos >= start_val,
                     pos <= stop_val)
  pheno_val <- pheno # need to be different from column name in x
  if(!is.null(pheno_val))
    x <-dplyr::filter(x, pheno %in% pheno_val)
  attr(x, "haplos") <- haplos
  class(x) <- unique(c("top_snps_tbl", class(x)))
  x
}

snpinfo_to_haplos <- function(snpinfo) {
  alleles <- dplyr::select(
    snpinfo,
    -(snp_id:alleles))
  infonames <- c("consequence","type","sdp","index","interval","on_map","pheno","lod")
  if(length(wh <- which(infonames %in% names(alleles)))) {
    alleles <- alleles[, -wh, drop = FALSE]
  }
  # Columns in between consequence and type should be alleles.
  # If not provided, assume we are in mouse with 8.
  if((nc <- ncol(alleles)) < 2) {
    warning("no alleles in snpinfo; assuming 8")
    nc <- 8
  }
  LETTERS[seq_len(nc)]
}

