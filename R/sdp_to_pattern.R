#' Convert sdp to pattern
#'
#' Convert strain distribution pattern (sdp) to letter pattern.
#'
#' @param sdp vector of sdp values
#' @param haplos letter codes for haplotypes (required)
#' @param symmetric make patterns symmetric if \code{TRUE}
#'
#' @return vector of letter patterns
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"
#' 
#' \donttest{
#' # Download SNP info for DOex from web and read as RDS.
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "c2_snpinfo.rds"), tmpfile, quiet=TRUE)
#' snpinfo <- readRDS(tmpfile)
#' unlink(tmpfile)
#' snpinfo <- dplyr::rename(snpinfo, pos = pos_Mbp)
#' 
#' # Extract strain distribution pattern.
#' sdp <- snpinfo$sdp
#' # Find out how many alleles.
#' nallele <- ceiling(log2(max(sdp)))
#' out <- sdp_to_pattern(sdp, LETTERS[seq_len(nallele)])
#' # Show most frequent patterns. 
#' head(rev(sort(c(table(out)))))
#' }
#'
#' @importFrom assertthat assert_that
#' @export
#' 
sdp_to_pattern <- function(sdp, haplos, symmetric = TRUE) {
  assertthat::assert_that(!missing(haplos))
  if(symmetric) {
    numsdp <- 2 ^ length(haplos)
    sdp <- ifelse(sdp < numsdp / 2, sdp, numsdp - sdp - 1)
  }
  sapply(sdp, function(x, haplos) {
    ref <- as.logical(intToBits(x)[seq_along(haplos)])
    paste(paste(haplos[!ref], collapse = ""),
          paste(haplos[ref], collapse = ""),
          sep = ":")
  }, haplos)
}
#' @rdname sdp_to_pattern
#' @export
sdp_to_logical <- function(sdp, haplos, symmetric = TRUE) {
  assertthat::assert_that(!missing(haplos))
  if(symmetric) {
    numsdp <- 2 ^ length(haplos)
    sdp <- ifelse(sdp < numsdp / 2, sdp, numsdp - sdp - 1)
  }
  sapply(sdp, function(x, haplos) {
    as.logical(intToBits(x)[seq_along(haplos)])
  }, haplos)
}
