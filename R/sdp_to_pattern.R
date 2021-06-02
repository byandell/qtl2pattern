#' Convert sdp to pattern
#'
#' Convert strain distribution pattern (sdp) to letter pattern.
#'
#' @param sdp vector of sdp values
#' @param haplos letter codes for haplotypes (required)
#'
#' @return vector of letter patterns
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' example(DOex)
#' # Extract strain distribution pattern.
#' sdp <- snpinfo$sdp
#' # Find out how many alleles.
#' nallele <- ceiling(log2(max(sdp)))
#' out <- sdp_to_pattern(sdp, LETTERS[seq_len(nallele)])
#' # Show most frequent patterns. 
#' head(rev(sort(table(out))))
#'
#' @export
#' @importFrom assertthat assert_that
#' 
sdp_to_pattern <- function(sdp, haplos) {
  assertthat::assert_that(!missing(haplos))
  sapply(sdp, function(x, haplos) {
    ref <- as.logical(intToBits(x)[seq_along(haplos)])
    paste(paste(haplos[!ref], collapse = ""),
          paste(haplos[ref], collapse = ""),
          sep = ":")
  }, haplos)
}
#' @export
#' @rdname sdp_to_pattern
sdp_to_logical <- function(sdp, haplos) {
  assertthat::assert_that(!missing(haplos))
  sapply(sdp, function(x, haplos) {
    as.logical(intToBits(x)[seq_along(haplos)])
  }, haplos)
}
