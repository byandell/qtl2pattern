#' Top SNPs organized by allele pattern
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
#' @param haplos optional argument identify codes for haplotypes
#'
#' @return table of top_snps at maximum lod for \code{pattern}
#'
#' @examples
#' dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"
#' 
#' # Read DOex example cross from 'qtl2data'
#' DOex <- subset(qtl2::read_cross2(file.path(dirpath, "DOex.zip")), chr = "2")
#' 
#' \donttest{
#' # Download genotype probabilities
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "DOex_genoprobs_2.rds"), tmpfile, quiet=TRUE)
#' pr <- readRDS(tmpfile)
#' unlink(tmpfile)
#' 
#' # Download SNP info for DOex from web and read as RDS.
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "c2_snpinfo.rds"), tmpfile, quiet=TRUE)
#' snpinfo <- readRDS(tmpfile)
#' unlink(tmpfile)
#' snpinfo <- dplyr::rename(snpinfo, pos = pos_Mbp)
#' 
#' # Convert to SNP probabilities
#' snpinfo <- qtl2::index_snps(DOex$pmap, snpinfo)
#' snppr <- qtl2::genoprob_to_snpprob(pr, snpinfo)
#' 
#' # Scan SNPs.
#' scan_snppr <- qtl2::scan1(snppr, DOex$pheno)
#' 
#' # Collect top SNPs
#' top_snps_tbl <- top_snps_pattern(scan_snppr, snpinfo)
#' summary(top_snps_tbl)
#' }
#' 
#' @export
#' @importFrom dplyr everything filter group_by inner_join select ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#'
top_snps_pattern <- function (scan1_output, snpinfo, drop = 1.5, show_all_snps = TRUE,
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
    lod_df <- tidyr::pivot_longer(
      dplyr::select(lod_df, .data$snp_id, .data$index, dplyr::everything()),
      -(1:2), names_to = "pheno", values_to = "lod")
    maxlod <- max(lod_df$lod, na.rm = TRUE) - drop
    lod_df <- dplyr::ungroup(
      dplyr::filter(
        dplyr::group_by(lod_df, .data$pheno),
        .data$lod >= min(max(.data$lod, na.rm = TRUE), maxlod)
      )
    )
#    lod_df <- dplyr::filter(lod_df, lod > max(lod, na.rm = TRUE) - drop)

    if (show_all_snps) {
      snpinfo <- dplyr::inner_join(snpinfo,
                                   dplyr::select(lod_df, -.data$snp_id),
                                   by = "index")
    }
    else {
      snpinfo <- dplyr::inner_join(snpinfo,
                                   dplyr::select(lod_df, -.data$index),
                                   by = "snp_id")
    }
    attr(snpinfo, "haplos") <- haplos
    class(snpinfo) <- c("top_snps_pattern", class(snpinfo))
    snpinfo
}

#' Summary of top SNP object across phenotypes
#'
#' @param object object of class \code{top_snps_tbl}
#' @param sum_type type of summary (one of "range","best")
#' @param ... other arguments not used
#'
#' @return table summary
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method summary top_snps_pattern
#' @rdname top_snps_pattern
#' @export
#' @importFrom dplyr arrange desc everything filter group_by mutate select summarize ungroup
#'
summary.top_snps_pattern <- function(object, sum_type=c("range","best","peak"), # peak is now old, not used
                                 ...) {
  haplos <- attr(object, "haplos")
  sum_type <- match.arg(sum_type)
  switch(sum_type,
         best = { ## Top SNPs across all phenotypes.
           out <- 
             dplyr::arrange(
               dplyr::select(object, -.data$index),
               dplyr::desc(.data$lod))
           if(!is.null(out$consequence)) {
             out <- 
               dplyr::mutate(out,
                 consequence = abbreviate(.data$consequence, 15))
           }
           dplyr::select(
             dplyr::group_by(
               dplyr::filter(
                 dplyr::group_by(out, .data$pheno, .data$sdp),
                 .data$lod == max(.data$lod))),
               .data$pheno,
               .data$chr, .data$pos,
               .data$lod,
               .data$snp_id, .data$sdp,
               dplyr::everything())
         },
         peak =,
         range = {
           ## SNP patterns within drop_hilit of max LOD.
           dplyr::select(
             dplyr::arrange(
               dplyr::mutate(
                 dplyr::ungroup(
                   dplyr::summarize(
                     dplyr::group_by(object, .data$sdp, .data$pheno),
                     max_n = ifelse(sum_type == "range",
                                    dplyr::n(),
                                    sum(.data$lod == max(.data$lod))),
                     min_pos = min(.data$pos[which(.data$lod == max(.data$lod))]),
                     max_pos = max(.data$pos[which(.data$lod == max(.data$lod))]),
                     snp_id = ifelse(.data$max_n == 1,
                                  .data$snp_id[which(.data$lod == max(.data$lod))][1],
                                  paste(.data$max_n, "SNPs")),
                     max_lod = max(.data$lod),
                     min_lod = min(.data$lod))),
                 pattern = sdp_to_pattern(.data$sdp, haplos)),
               dplyr::desc(.data$max_lod)),
             .data$pheno, 
             .data$min_pos, .data$max_pos,
             .data$max_lod, .data$min_lod,
             .data$sdp, .data$pattern, .data$snp_id)
         })
}
#' Subset of features
#'
#' @param x tbl of feature information from \code{\link{get_feature_snp}}
#' @param start_val,end_val start and end positions for subset
#' @param pheno phenotype name(s) for subset
#' @param ... additional parameters ignored
#'
#' @return subset of \code{x}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method subset top_snps_pattern
#' @rdname top_snps_pattern
#' @export
#' @importFrom dplyr filter
subset.top_snps_pattern <- function(x, start_val=0, end_val=max(x$pos),
                                pheno = NULL, ...) {
  haplos <- attr(x, "haplos")
  x <- dplyr::filter(x,
                     .data$pos >= start_val,
                     .data$pos <= end_val)
  pheno_val <- pheno # need to be different from column name in x
  if(!is.null(pheno_val))
    x <-dplyr::filter(x, .data$pheno %in% pheno_val)
  attr(x, "haplos") <- haplos
  class(x) <- unique(c("top_snps_tbl", class(x)))
  x
}
